import cobra
from cobra.io import load_json_model, save_json_model
from cobra import Model, Reaction, Metabolite
import numpy as np
from math import exp
import pandas as pd


ex_aa_list = ['EX_asn__L_e', 'EX_cys__L_e', 'EX_gln__L_e', 'EX_lys__L_e', 'EX_pro__L_e', 
              'EX_tyr__L_e', 'EX_met__L_e', 'EX_leu__L_e', 'EX_ser__L_e', 'EX_his__L_e', 
              'EX_thr__L_e', 'EX_phe__L_e', 'EX_arg__L_e', 'EX_ile__L_e', 'EX_val__L_e', 
              'EX_trp__L_e', 'EX_asp__L_e', 'EX_ala__L_e', 'EX_glu__L_e', 'EX_gly_e']

def set_PA(model, ptot, A_dict):
    '''
    Set proteome constraints to a GEM.
    model: GEM model
    ptot: total protein mass per 1 gram dry weight biomass
    A_dict: the dict of enzyme activities
    '''
    sigma=0.5
    t_sector = model.reactions.EX_lac__L_e.flux_expression/(sigma*A_dict['EX_lac__L_e']) +\
           (-1)*model.reactions.EX_glc__D_e.flux_expression/A_dict['EX_glc__D_e'] +\
             model.reactions.EX_ac_e.flux_expression/(sigma*A_dict['EX_ac_e'])
    for ex_aa in ex_aa_list:
        if ex_aa != 'EX_ala__L_e':
            t_sector = t_sector + (-1)*model.reactions.get_by_id(ex_aa).flux_expression/A_dict['EX_aa_e']
    a_sector = model.reactions.BIOMASS.flux_expression/A_dict['BIOMASS']
    ngam_sector = model.reactions.ATPM.flux_expression/A_dict['ATPM']
    c_sector = model.reactions.ENO.flux_expression/(sigma*A_dict['ENO'])
    for k in A_dict.keys():
        if k not in ['ENO','EX_lac__L_e', 'EX_ac_e','EX_glc__D_e','EX_aa_e','ATPM','BIOMASS']:
            c_sector  = c_sector  + model.reactions.get_by_id(k).flux_expression/(sigma*A_dict[k])
    PA = model.problem.Constraint( expression = a_sector + c_sector + ngam_sector,
                        name = 'PA', lb= 0, ub = ptot)
    model.add_cons_vars([ PA ])
    return None

def pcfba(model, ptot, NGAM, AA_lb, Glc_lb, OUR_lb, A_dict):
    '''
    Run proteome constrained flux balance analysis (for PigGEM2025).
    model: GEM model
    ptot: total protein mass per 1 gram dry weight biomass
    NGAM: non-growth associated maintenance ATP
    AA_lb: (negative) lower bound of amino acid uptake
    Glc_lb: (negative) lower bound of glucose uptake
    OUR_lb: (negative) lower bound of oxygen uptake
    A_dict: the dict of enzyme activities
    '''
    with model:
        for ex_aa in ex_aa_list:
            model.reactions.get_by_id(ex_aa).lower_bound=-AA_lb
        model.reactions.EX_glc__D_e.lower_bound=-Glc_lb
        model.reactions.ATPM.lower_bound = NGAM
        model.reactions.EX_o2_e.lower_bound = -OUR_lb
        set_PA(model, ptot, A_dict)
        fluxes = cobra.flux_analysis.pfba(model).fluxes
    return fluxes

def lac_inhibit_glc(lac_con, v_min):
    #vmin=0.01
    return max( -0.084*lac_con + 1.789, v_min )

def nh4_inhibit_glc(lac_con, nh4_con, v_min):
    #vmin=0.01
    if nh4_con <= 4:
        return lac_inhibit_glc(lac_con, v_min)
    return max( lac_inhibit_glc(lac_con, v_min)*0.9612/(nh4_con-2.9434), v_min )

def inhibit_gln(lac_con, nh4_con, v_min):
    #vmin=0.01
    k1,k2,k3,k4 = 4.688,2.346,5.544,2.907
    return max( 0.526*(k1/(lac_con+k2))*(k3/(nh4_con+k4)), v_min )

def inhibit_our(lac_con, nh4_con, v_min):
    #vmin=0.6
    k1,k2,k3,k4 = 10.265,1.512,1.812,2.255
    return max( 4.9*(k1/(lac_con+k2))*(k3/(nh4_con+k4)), v_min )

def dpcfba(model, ptot, NGAM, Dr, A_dict, ic, T):
    '''
    Dynamic proteome constrained flux balance analysis
    model: GEM model
    ptot: total protein mass per 1 gram dry weight biomass
    NGAM: non-growth associated maintenance ATP
    Dr: death rate, 0.01/hr
    A_dict: the dict of enzyme activities
    ic: initial condition
    T: time length (hr)
    '''
    times,step = np.linspace(0,T,num= 100,retstep=True)
    met_profile = {key: [ic[key]] for key in ic.keys() }
    flux_profile = {'BIOMASS':[]}
    for k in met_profile:
        if k != 'BIOMASS':
            flux_profile['EX_'+k] = []
            
    for i in range(len(times)-1):
        profile_t = {k: met_profile[k][i] for k in met_profile.keys() }#concentration profile snapshot at time t = i
        AA_lb = inhibit_gln(profile_t['lac__L_e'], profile_t['nh4_e'], 0.01)
        Glc_lb = nh4_inhibit_glc(profile_t['lac__L_e'], profile_t['nh4_e'], 0.01)
        OUR_lb = inhibit_our(profile_t['lac__L_e'], profile_t['nh4_e'], 0.6)
        fluxes = pcfba(model, ptot, NGAM, AA_lb, Glc_lb, OUR_lb, A_dict)
        for k in flux_profile:
            flux_profile[k].append( fluxes[k] )
        for k in met_profile.keys():
            if k == 'BIOMASS':
                met_profile['BIOMASS'].append( met_profile['BIOMASS'][i]+\
                                              met_profile['BIOMASS'][i]*(flux_profile['BIOMASS'][i]-Dr)*step)
            else:
                met_profile[k].append( max(met_profile[k][i]+met_profile['BIOMASS'][i]*flux_profile['EX_'+k][i]*step,0) )#conc >= 0
    met_profile['T'] = times
    return met_profile, flux_profile
                
    
    
    

    
            
        
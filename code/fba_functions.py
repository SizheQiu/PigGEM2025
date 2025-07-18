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
    t_sector = model.reactions.EX_lac__L_e.flux_expression/A_dict['EX_lac__L_e'] +\
           (-1)*model.reactions.EX_glc__D_e.flux_expression/A_dict['EX_glc__D_e'] +\
             model.reactions.EX_ac_e.flux_expression/A_dict['EX_ac_e']
    for ex_aa in ex_aa_list:
        t_sector = t_sector + (-1)*model.reactions.get_by_id(ex_aa).flux_expression/A_dict['EX_aa_e']
    a_sector = model.reactions.BIOMASS.flux_expression/A_dict['BIOMASS']
    ngam_sector = model.reactions.ATPM.flux_expression/A_dict['ATPM']
    c_sector = model.reactions.ENO.flux_expression/A_dict['ENO']
    for k in A_dict.keys():
        if k not in ['ENO','EX_lac__L_e', 'EX_ac_e','EX_glc__D_e','EX_aa_e','ATPM','BIOMASS']:
            c_sector  = c_sector  + model.reactions.get_by_id(k).flux_expression/A_dict[k]
    PA = model.problem.Constraint( expression = a_sector + c_sector + ngam_sector,
                        name = 'PA', lb= 0, ub = ptot)
    model.add_cons_vars([ PA ])
    return None

def pcfba(model, ptot, NGAM, AA_lb, Glc_lb, A_dict):
    with model:
        for ex_aa in ex_aa_list:
            model.reactions.get_by_id(ex_aa).lower_bound=-AA_lb
        model.reactions.EX_glc__D_e.lower_bound=-Glc_lb
        model.reactions.ATPM.lower_bound = NGAM
        set_PA(model, ptot, A_dict)
        fluxes = cobra.flux_analysis.pfba(ipsc_core).fluxes
    return fluxes

def lac_inhibit_glc(lac_con, v_min):
    return min( -0.084*lac_con + 1.789, v_min )

def nh4_inhibit_glc(lac_con, nh4_con, v_min):
    if nh4_con <= 4:
        return lac_inhibit_glc(lac_con, v_min)
    return min( lac_inhibit_glc(lac_con, v_min)*2.299/(nh4_con-0.146), v_min )

    
            
        
#!/bin/bash
#SBATCH --ntasks-per-node=1
#SBATCH --job-name=blast
#SBATCH --partition=short
#SBATCH --time=10:00:00
#SBATCH --mem-per-cpu=32G


module load BLAST+/2.7.1-foss-2018b

blastp -query recon1_protein.fasta -subject ss_protein.faa -evalue 1e-6 -outfmt 7 -out human_blast_ss.txt
blastp -query icho_protein.fasta -subject ss_protein.faa -evalue 1e-6 -outfmt 7 -out cho_blast_ss.txt

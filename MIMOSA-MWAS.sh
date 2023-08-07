#!/bin/bash                                                                    
#SBATCH -J MWAS
#SBATCH -n 2
#SBATCH -p backfill2
#SBATCH -t 4:00:00
#SBATCH --array=1-300
#SBATCH --mail-type=ALL
#SBATCH --mail-user=YOUR EMAIL HERE

#SBATCH -C "YEAR2017|YEAR2018|YEAR2019"                                         

export PATH=$PATH:/gpfs/research/chongwu/shared/software

module load gnu/9.1.1
module load gnu-openmpi

module load R/4.2.0

R CMD BATCH --no-save --no-restore "--args path.ref='PATH TO LD REF' trait='TRAIT ABBREVIATION/NAME' path.trait='PATH TO GWAS SUMMARY STATS' path.out='PATH TO SAVE RESULTS' path.weight='PATH TO WHERE MIMOSA MODELS ARE SAVED'" MIMOSA-MWAS.R ./log/MWAS_$SLURM_ARRAY_TASK_ID.txt


#!/bin/bash
#SBATCH --job-name=run_hmm
#SBATCH --time=05:00:00
#SBATCH --ntasks=1
#SBATCH -p hns,schumer,owners
#SBATCH --cpus-per-task=1
#SBATCH --mem=32000
#SBATCH --mail-user=schumer@stanford.edu

module load python/3.9
module load boost/1.76.0
module load armadillo
module load biology
module load samtools
module load bcftools
module load bwa
module load gcc/10.1.0
module load gsl
ml R
ml java
export PATH="/home/groups/schumer/shared_bin/Ancestry_HMM/src:$PATH"
export PATH="/home/groups/schumer/shared_bin:$PATH"

perl /home/groups/schumer/shared_bin/ancestryinfer_July2022/Ancestry_HMM_parallel_v7.pl hmm_configuration_file.cfg

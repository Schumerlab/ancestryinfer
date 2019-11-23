#!/bin/bash
#SBATCH --job-name=run_hmm
#SBATCH --time=05:00:00
#SBATCH --ntasks=1
#SBATCH -p schumer
#SBATCH --cpus-per-task=1
#SBATCH --mem=32000
#SBATCH --mail-user=schumer@stanford.edu

module load armadillo
module load biology
module load samtools
module load bcftools
module load py-pysam/0.14.1_py27
module load bwa
export PATH="/home/groups/schumer/shared_bin/ngsutils/bin:$PATH"
export PATH="/home/groups/schumer/shared_bin/Ancestry_HMM/src:$PATH"
export PYTHONPATH=/home/groups/schumer/shared_bin:$PYTHONPATH

perl  /home/groups/schumer/shared_bin/Ancestry_HMM_pipeline/Ancestry_HMM_parallel_v5.pl hmm_configuration_file.cfg


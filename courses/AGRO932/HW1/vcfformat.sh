#!/bin/bash -l
#SBATCH --output=/common/soybean/bharms4/slurm-log/output/vcfHW1soy-stdout-%j.txt
#SBATCH --error=/common/soybean/bharms4/slurm-log/error/vcfHW1soy-stderr-%j.txt
#SBATCH --job-name=vcfHW1soy
#SBATCH --time 24:00:00
#SBATCH --mail-user=benjamin.harms@huskers.unl.edu
#SBATCH --mail-type=ALL #email if ends
#SBATCH --mail-type=FAIL #email if fails

module load bwa samtools
cd /common/soybean/bharms4/courses/2022-agro932-lab/largedata1/HW1

samtools mpileup -g -f glycinemaxMT.fa -b bamlist.txt > myraw.bcf


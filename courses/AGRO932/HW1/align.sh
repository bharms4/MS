#!/bin/bash -l
#SBATCH --output=/common/soybean/bharms4/slurm-log/output/alignHW1soy-stdout-%j.txt
#SBATCH --error=/common/soybean/bharms4/slurm-log/error/alignHW1soy-stderr-%j.txt
#SBATCH --job-name=alignHW1soy
#SBATCH --time 24:00:00
#SBATCH --mail-user=benjamin.harms@huskers.unl.edu
#SBATCH --mail-type=ALL #email if ends
#SBATCH --mail-type=FAIL #email if fails

module load bwa samtools
cd /common/soybean/bharms4/courses/2022-agro932-lab/largedata1/HW1
# alignment
for i in {1..20}; do bwa mem glycinemaxMT.fa l$i.read1.fq l$i.read2.fq | samtools view -bSh - > l$i.bam; done
# sort
for i in *.bam; do samtools sort $i -o sorted_$i; done
# index them
for i in sorted*.bam; do samtools index $i; done

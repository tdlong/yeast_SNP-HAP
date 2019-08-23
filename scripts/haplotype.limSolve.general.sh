#!/bin/bash
#$ -N lhap
#$ -q bio,pub8i,abio128
#$ -tc 600
#$ -ckpt restart
#$ -pe openmp 2
#$ -R y

module load R/3.4.1
file=$1
chr=`head -n $SGE_TASK_ID $file | tail -n 1 | cut -f1` 
pool=`head -n $SGE_TASK_ID $file | tail -n 1 | cut -f2` 
# haplotyper.4.code.R is current the 50:50 caller with sigma
# chosen so 50 sites closest to position
# account for 50% of weight
Rscript scripts/haplotyper.limSolve.R $3 $chr $pool $4 $2


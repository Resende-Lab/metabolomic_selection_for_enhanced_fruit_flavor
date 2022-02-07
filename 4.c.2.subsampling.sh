#!/bin/bash
#SBATCH --job-name=9.tmp/2.subsamp.%j
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --nodes=1
#SBATCH --mem=4GB
#SBATCH -t 48:00:00
#SBATCH --output=9.tmp/2.subsamp.%a.%A.out
#SBATCH --array=1-1

module load R
module load ufrc

experiment=$1
trait=$2
subsamp=$3

echo $experiment $trait $subsamp
Rscript 4.c.3.subsampling.R $experiment $trait $subsamp

#!/bin/bash
#SBATCH --job-name=9.tmp/3.jobs/array.%j
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --nodes=1
#SBATCH --mem=4GB
#SBATCH --qos=mresende-b
#SBATCH --account=mresende
#SBATCH -t 48:00:00
#SBATCH --output=9.tmp/3.jobs/%a.array.%A.out
#SBATCH --array=1-1

module load R

experiment=$1
trait=$2
cv=$3

echo $experiment $trait $cv
Rscript 4.b.6.metabolomic_selection_tomato.R $experiment $trait $cv


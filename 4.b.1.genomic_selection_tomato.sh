#!/bin/bash
#SBATCH --job-name=9.tmp/3.jobs/array.%j
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --nodes=1
#SBATCH --mem=4GB
#SBATCH --qos=mresende-b
#SBATCH --account=mresende
#SBATCH -t 1:00:00
#SBATCH --output=9.tmp/3.jobs/%a.array.%A.out
#SBATCH --array=1-1


for experiment in {1..100}
  do
  	for trait in {1..5}
  	do
  		for cv in {1..10}
  		do
  			echo $experiment $trait $cv
  			sbatch 4.b.2.genomic_selection_tomato.sh $experiment $trait $cv
  		done
  	done
done


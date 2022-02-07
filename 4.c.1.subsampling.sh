#!/bin/bash
#SBATCH --job-name=9.tmp/1.subsamp.%j
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --nodes=1
#SBATCH --mem=4GB
#SBATCH -t 1:00:00
#SBATCH --output=9.tmp/1.subsamp.%a.%A.out
#SBATCH --array=1-1

module load ufrc

for experiment in {1..10}
  do
  	for trait in {1..5}
  	do
  		for subsamp in {50..170..10}
  		do
  			echo $experiment $trait $subsamp
  			sbatch 4.c.2.subsampling.sh $experiment $trait $subsamp
  		done
  	done
done


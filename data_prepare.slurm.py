# -*- coding: utf-8 -*-
"""
Created on Fri Jan 22 12:08:31 2021

@author: Max Jansen
"""

#!/bin/bash
#SBATCH --nodes 1
#SBATCH --partition=serial
#SBATCH --ntasks-per-node 1
#SBATCH --cpus-per-task 1
#SBATCH --mem 8192
#SBATCH --time 01:00:00
#SBATCH --array=20-1000
#SBATCH --output=exelogs/interface_hydrophobicity.%A_%a.out
#SBATCH --error=exelogs/interface_hydrophobicity.%A_%a.err

i=1
while read p; do
    if [ $(( i % 1001 )) == ${SLURM_ARRAY_TASK_ID} ]; then
        FIELD1=$(echo $p| cut -d" " -f1)
        ./interface_calculator.py $FIELD1
    fi
    i=$((i+1))
done < list_ab_train.dat

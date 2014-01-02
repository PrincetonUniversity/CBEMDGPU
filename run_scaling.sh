#!/bin/bash
#PBS -l nodes=1:ppn=4,walltime=0:10:00

module load intel

cd $PBS_O_WORKDIR

for i in 1 2 4
  do
  for j in 200 400 800 1600
    do
    for k in 0 0.1 0.2 
      do
      ./scaling $i $j $k 10000
    done
  done
done
#!/bin/csh

#PBS -l walltime=99:00:00
#PBS -l nodes=1:ppn=1
#PBS -j oe


module add intel-compilers/12.0.4.191

cd $PBS_O_WORKDIR
./main



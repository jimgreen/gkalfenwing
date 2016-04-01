#!/bin/bash
#PBS -N m3a01
#PBS -l nodes=6:ppn=32
#PBS -l walltime=24:00:00
#PBS -r n
#PBS -j oe

NPROCS=`wc -l < $PBS_NODEFILE`
NNODES=`uniq $PBS_NODEFILE | wc -l`

EXE="/home/space/liang/software/gkeyll-par/bin/gkeyll"
FILE="mirdip.lua"

cd $PBS_O_WORKDIR
aprun -n $NPROCS $EXE -i $FILE  2>&1 | tee log



#!/bin/bash
echo "Bash version ${BASH_VERSION}"
mpicc diff_expl_fd.c -o diffusion -lm
echo "Enter grid points along x:"
read gridx
echo "Enter grid points along y:"
read gridy

declare -a TIME 
count=0

for i in 1 2 4 8 16 32
do
procs=$i
echo "Running mpiexec -n $procs ./diffusion $gridx $gridy"
STARTTIME=$(date +%s)
mpiexec -n $procs ./diffusion $gridx $gridy
ENDTIME=$(date +%s)
TIME[$count]=$(($ENDTIME-$STARTTIME))
count=$(($count+1))
done

count=0
for i in 1 2 4 8 16 32
do
procs=$i
echo "Elapsed ${TIME[$count]} seconds with $procs processor"
count=$(($count+1))
done

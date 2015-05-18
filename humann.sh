#!/bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -V
#$ -M wh359@drexel.edu
#$ -P nsftuesPrj
#$ -l h_rt=4:00:00
#$ -l h_vmem=8G
#$ -q all.q@@amdhosts
#$ -pe shm 2-16

 /etc/profile.d/modules.sh
module load shared
module load proteus
module load sge/univa
module load gcc/4.8.1
module load ncbi-blast/gcc/64/2.2.29

cd /home/wh359/humann-0.99

scons

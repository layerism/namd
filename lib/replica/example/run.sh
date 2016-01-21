#!/bin/sh

namd_replica=/home/luyukun/software/app/namd/NAMD_CVS-2013-10-17_Source/Linux-x86-g++.arch/namd2

rm -rf output
mkdir output
cd output
rm -rf 0 1
mkdir 0 1

cd ../

mpirun -np 4 namd2  pair.namd > pair.log
#mpirun -np 4  $namd_replica +replicas 2 job0.conf +stdout output/%d/job0.%d.log
#mpirun -np 4  $namd_replica +replicas 4 job0.conf  output/%d/job0.%d.log
mpirun -np 4  $namd_replica  alanin_base.namd > md.log


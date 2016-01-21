#!/usr/bin/sh



rm -rf Linux-x86-g++.arch
./config tcl fftw Linux-x86-g++.arch --charm-base /PATH/TO/YOUR/OWN/NAMD/NAMD_CVS-2013-10-17_Source/charm-6.5.1 \
--charm-arch mpi-linux-gcc


cd Linux-x86-g++.arch
make
cd ../

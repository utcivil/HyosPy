#! /bin/bash


mpirun -np 2 --hostfile machines pelfe_3.1dc_sirius
./autocombine_MPI_elfe.pl 169 216 0 2

#!/bin/bash

# e4shared --job=nonaka --run --report-residuals  --verbosity=1 | tee LOGFILE 
mpirun -np 16 e4mpi --job=nonaka --run --report-residuals  --verbosity=1 | tee LOGFILE 

# e4shared --job=nonaka --tindx-start=5 --run --report-residuals \
#  --verbosity=1 | tee -a LOGFILE
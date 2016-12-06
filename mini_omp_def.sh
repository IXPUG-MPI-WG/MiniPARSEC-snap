#!/bin/bash

NTHREADS=${1:-6}
SCHED=${2:-"guided,1000"}
AFFINITY=${3:-"close"}

export OMP_STACKSIZE=500M #should work now
export OMP_SCHEDULE=$SCHED
export OMP_NUM_THREADS=$NTHREADS
export OMP_PROC_BIND=$AFFINITY
export OMP_PLACES=cores
export OMP_WAIT_POLICY=active
export OMP_DYNAMIC=.FALSE.

export MKL_NUM_THREADS=$NTHREADS
export MKL_DYNAMIC=false

export KMP_SETTINGS=1

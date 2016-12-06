miniPARSEC
==========

miniPARSEC is an MPI3 finite difference benchmark that mimics the halo exchange
(ghost exchange) implemented in PARSEC. 

It now has a hybrid MPI3-openMP4.0 version that is recommended for all runs.
See instructions below on how to run it 'properly'

Dev status - perpetual alpha stage. 

##Introduction

TBD, something about this being example of use of MPI 3.0 with openmp etc.


##Building

Edit UserSettings.mk , and make. 

The recommendation is to use `make all -jXX` or `make non-itac -jXX` for building everything or everything not including the ITAC libs.
(Good compile times were achieved when XX=8)

##How to run this benchmark

The miniparsec binary takes either 0, 2 or 4 options. No input options 
means that the defaults will be used. 

The bash script mini_def_run.sh shows an example of running all the binaries compiled
with all the options, and includes an option to run the ITAC binaries as well (if found). 
For A quick analysis, try `grep loop: */miniparsec.out`

If running the openMP version, please source the mini_omp_def.sh file first,
an optional argument (default 6) is the number of threads per mpi process.

Note, that it is up to the user to make sure that intel mpi pins the process to the correct NUMA node.
I suggest using I_MPI_PIN_DOMAIN=omp, which deals with it for you such as long $n_{mpi}*n_{omp}=n_{cpu}$.
Note that this may not be supported on certain clusters lacking memory affinity modules.

##Contributors

 1. Ariel Biller, Weizmann Institute of Science
 2. Charles Lena, University of Texas
 3. *Add your name here*

##License

GPL vX, see LICENSE file @todo

##Contact

ariel.biller@weizmann.ac.il

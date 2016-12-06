# fortran compiler
F90=mpiifort

# compiler openmp flag
OPENMP_FLAG=-qopenmp

# flag for directing where modules go to live (used as $(MODULE_PATH) $(PATH))
MODULE_PATH=-module

# for building ITAC examples - others will build without it
VT_INCLUDES=-I$(VT_ROOT)/include/
VT_LIB_LINK=-L$(VT_LIB_DIR) -lVT $(VT_ADD_LIBS) 
LIBMPI=$(VT_LIB_LINK)
# needs dlarnv_ for initalization 
LAPACK_LIB=-mkl

# base compiler optimization + VT support
FCFLAGS=-O3 -align array64byte -trace -traceback -ip $(VT_INCLUDES)

TYPE_PROFILER=itac
TYPE_BC=nbc
TYPE_OMP=omp
CPPOPT+=-DITAC -DOMP
EXT=.mpi

ifeq ($(origin VT_INCLUDES),undefined)
$(error VT_INCLUDES "must be defined in your environment or in UserSettings.mk to compile ITAC examples")
endif

FCFLAGS+=$(OPENMP_FLAG) $(PROFILER_FLAGS) $(VT_INCLUDES)

ifeq ($(origin VT_LIB_LINK),undefined)
$(error VT_LIB_LINK "must be define in your environment or in UserSettings.mk to link ITAC required libraries e.g. VT_LIB_LINK=-L$$(VT_ROOT) -lVT <morelibs here>.")
endif

LIBS+=$(VT_LIB_LINK)

# Source
SRC=mini_benchmark_mod.F90 mini_buffers_mod.F90 mini_common_mod.F90 mini_const_mod.F90 mini_global_data_mod.F90 mini_grid_mod.F90 mini_grid_partition.F90 mini_io_mod.F90 mini_parallel_data_mod.F90 mini_setup.F90 miniparsec.F90
OBJ=$(SRC:.F90=.o)

# Takes the configuration and makes it into a name
TYPES=$(if $(strip $(TYPE_PROFILER)),.,)$(strip $(TYPE_PROFILER))$(if $(strip $(TYPE_BC)),.,)$(strip $(TYPE_BC))$(if $(strip $(TYPE_OMP)),.,)$(strip $(TYPE_OMP))

# Represents the build directory
BUILDDIR=obj/$(TYPES)
TMPSRCDIR=tmp/$(TYPES)
BUILD_OBJ=$(addprefix $(BUILDDIR)/,$(POBJF) $(POBJZ) $(OBJ))
$(shell mkdir -p $(BUILDDIR))
$(shell mkdir -p $(TMPSRCDIR))

miniparsec$(TYPES)$(EXT): $(BUILD_OBJ)
	    $(F90) -o $@ $(FCFLAGS) $^ $(LAPACK_LIB) $(LIBMPI) $(LIBS)

$(BUILDDIR)/%.o: %.f
	 $(FC) -c $(FCFLAGS) $*.f $(MODULE_PATH) $(BUILDDIR) -o $(BUILDDIR)/$*.o

$(BUILDDIR)/%.o: %.f90
	 $(F90) -c $(FCFLAGS) $*.f90 $(MODULE_PATH) $(BUILDDIR) -o $(BUILDDIR)/$*.o

$(BUILDDIR)/%.o: %.F90
	 $(F90) $(CPPOPT) -c $(FCFLAGS) $(MODULE_PATH) $(BUILDDIR) $< -o $(BUILDDIR)/$*.o

$(BUILDDIR)/%.o: %.f90p
	 $(CPP) $(CPPOPT) $< | sed '/SOME_SED_CALL/d' > $(TMPSRCDIR)/$*.tmp.f90
	 $(F90) -c $(FCFLAGS) $(MODULE_PATH) $(TMPSRCDIR)/$*.tmp.f90 -o $(BUILDDIR)/$*.o

$(BUILDDIR)/%.o: %.f90z
	 $(CPP) $(CPPOPT) $< | sed '/SOME_SED_CALL/d' > $(TMPSRCDIR)/$*.tmp.f90
	 $(CPP) $(CPPOPT) -DCPLX $< | sed '/CPLX_SED_CALL/d' >> $(TMPSRCDIR)/$*.tmp.f90
	 $(F90) -c $(FCFLAGS) $(MODULE_PATH) $(TMPSRCDIR)/$*.tmp.f90 -o $(BUILDDIR)/$*.o

# SOURCE DEPENDENCIES FOR MODULES
$(BUILDDIR)/mini_parallel_data_mod.o $(BUILDDIR)/mini_common_mod.o: $(BUILDDIR)/mini_const_mod.o 
$(BUILDDIR)/mini_grid_mod.o: $(BUILDDIR)/mini_common_mod.o
$(BUILDDIR)/mini_global_data_mod.o: $(BUILDDIR)/mini_grid_mod.o $(BUILDDIR)/mini_parallel_data_mod.o $(BUILDDIR)/mini_buffers_mod.o
$(BUILDDIR)/mini_buffers_mod.o:  $(BUILDDIR)/mini_parallel_data_mod.o
$(BUILDDIR)/mini_grid_partition.o:  $(BUILDDIR)/mini_parallel_data_mod.o $(BUILDDIR)/mini_global_data_mod.o
$(BUILDDIR)/mini_setup.o: $(BUILDDIR)/mini_common_mod.o $(BUILDDIR)/mini_const_mod.o $(BUILDDIR)/mini_global_data_mod.o $(BUILDDIR)/mini_buffers_mod.o
$(BUILDDIR)/mini_io_mod.o: $(BUILDDIR)/mini_global_data_mod.o
$(BUILDDIR)/mini_benchmark_mod.o: $(BUILDDIR)/mini_common_mod.o $(BUILDDIR)/mini_const_mod.o $(BUILDDIR)/mini_global_data_mod.o $(BUILDDIR)/mini_buffers_mod.o
$(BUILDDIR)/miniparsec.o: $(BUILDDIR)/mini_benchmark_mod.o

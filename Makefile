#!/usr/bin/env make
# general compile choices go here
include UserSettings.mk

ifneq ($(strip $(TARGET)),)
# magic beans from giants
include Makefile.$(TARGET).mk
# source file changes live here
include BuildRules.mk
endif
default:
	@echo ""
	@echo "" choose:
	@echo ""
	@echo "" make nbc      = non-blocking neighborhood collectives
	@echo "" make nbc-omp  = non-blocking neighborhood collectives with basic openMP kernel
	@echo "" make mpi2     = MPI2-style irecv/isend+waitall
	@echo "" make omp-mpi2 = MPI2-style irecv/isend+waitall with basic openMP kernel
	@echo "" make bc       = blocking neighborhood collectives
	@echo "" make bc-omp   = blocking neighborhood collectives with basic openMP kernel
	@echo ""
	@echo "" make itac-nbc      += with ITAC api
	@echo "" make itac-nbc-omp  += with ITAC api
	@echo "" make itac-mpi2     += with ITAC api
	@echo "" make itac-mpi2-omp += with ITAC api
	@echo "" make itac-bc       += with ITAC api
	@echo "" make itac-bc-omp   += with ITAC api
	@echo ""
	@echo "" or use:
	@echo ""
	@echo "" make all      \(makes everything\)
	@echo "" make non-itac \(makes bc,nbc,mpi2 and their omp variants\)
	@echo "" make all-itac \(makes bc,nbc,mpi2 and their omp variants\) with ITAC
	@echo ""
	@echo ""
	@echo "" Tips:
	@echo ""
	@echo "" make clean/veryclean = object/all cleanup
	@echo "" make cleantest = clean directories created by the mini_def_run.sh script
	@echo "" make tags = generate ctags for this project
	@echo "" make doxy = generaye documentation using doxygen - beta
	@echo "" 

tags:
	ctags --if0=yes --language-force=Fortran -R *.F90 -h ".inc" 

clean:
	rm -f *.o *.mod
	rm -rf obj/ tmp/

veryclean: clean
	rm -f miniparsec*.mpi

doxy:
	doxygen Doxyfile.conf

.PHONY: clean veryclean tags nbc nbc-omp bc bc-omp mpi2 mpi2-omp itac-nbc itac-nbc-omp itac-bc miniparsec
.SUFFIXES: .f90p .f90z .F90 .f90

PROGRAMS=nbc \
	 nbc-omp \
	 bc \
	 bc-omp \
	 mpi2 \
	 mpi2-omp \
	 itac-nbc \
	 itac-nbc-omp \
	 itac-bc \
	 itac-bc-omp \
	 itac-mpi2 \
	 itac-mpi2-omp

all: $(PROGRAMS) fxtest
all-itac: itac-bc itac-bc-omp itac-nbc-omp itac-nbc itac-mpi2 itac-mpi2-omp fxtest
non-itac: nbc nbc-omp bc bc-omp mpi2 mpi2-omp fxtest


# SOURCE BUILD RULES
miniparsecs: 
	$(MAKE) miniparsec$(TYPES)$(EXT)
	@echo $(TARGET) finished!

# UTILITY RULES
fxtest:
	chmod +x mini_def_run.sh

cleantest:
	@echo "      cleaning directories created by mini_def_run.sh"
	@rm t_*/ trace-tmp/ -r 2>/dev/null

############### rule that requires a variable to be defined or else we exit
define require
ifneq ($(strip $(2)),)
$(1)_msg = $(2)
else
$(1)_msg = $(1) should really be defined. Missing message.
endif
$(error $(1) - $($(1)_msg))
endef

############# The following auto builds rules that look like this:
# nbc:  
# 	echo $@: $(TYPES)
# 	$(MAKE) TARGET=nbc miniparsecs
################
# INPUT: prog in PROGRAMS
# OUTPUT: Text rules for make
define prog_MAKE
$(1): Makefile.$(1).mk
	echo $$@: $$(TYPES)
	$$(MAKE) TARGET=$$@ miniparsecs

clean-$(1):
	echo $$@: $$(TYPES)
	rm -rf ./obj/$$(TYPES) ./tmp/$$(TYPES)

endef
$(foreach prog,$(PROGRAMS),$(eval $(call prog_MAKE,$(prog))))


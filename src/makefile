# This makefile uses intel fortran compiler for linux

# Makefile variables
FC      = ifort
FCFLAGS = -g -parallel -mkl -traceback -CB
VPATH	:= src
OBJDIR	:= ../obj
LIBPATH := /opt/intel/mkl/lib/intel64
LIBS    := $(LIBPATH)/libmkl_intel_lp64.a \
  			$(LIBPATH)/libmkl_core.a  \
  			$(LIBPATH)/libmkl_sequential.a \
  			$(LIBPATH)/libmkl_core.a  \
  			$(LIBPATH)/libmkl_sequential.a \
  			$(LIBPATH)/libmkl_core.a  \
			$(LIBPATH)/libmkl_intel_ilp64.a \
  			-lpthread

# Create the object directory to hold all .o, .mod and executables
$(OBJDIR):
	mkdir -p $(OBJDIR)

# All objects
OBJECTS := 	$(OBJDIR)/runModel.o $(OBJDIR)/fem.o \
		   	$(OBJDIR)/rt.o $(OBJDIR)/mesh.o \
		   	$(OBJDIR)/femHelper.o $(OBJDIR)/rtHelper.o \
			$(OBJDIR)/commonRoutines.o

# Make all
all: $(OBJECTS)
	$(FC) $(FCFLAGS) -o $(OBJDIR)/test $(OBJECTS) $(LIBS)

$(OBJDIR)/runModel.o: 	$(OBJDIR)/fem.o $(OBJDIR)/rt.o \
						$(OBJDIR)/mesh.o $(OBJDIR)/femHelper.o \
						$(OBJDIR)/rtHelper.o $(OBJDIR)/commonRoutines.o | $(OBJDIR)
	$(FC) $(FCFLAGS) -I$(OBJDIR) -c runModel.f90 -o $(OBJDIR)/runModel.o

$(OBJDIR)/fem.o: 	$(OBJDIR)/mesh.o $(OBJDIR)/femHelper.o \
					$(OBJDIR)/commonRoutines.o | $(OBJDIR)
	$(FC) $(FCFLAGS) -I$(OBJDIR) -module $(OBJDIR) -c fem.f90 -o $(OBJDIR)/fem.o

$(OBJDIR)/rt.o: 	$(OBJDIR)/mesh.o $(OBJDIR)/rtHelper.o \
					$(OBJDIR)/commonRoutines.o | $(OBJDIR)
	$(FC) $(FCFLAGS) -I$(OBJDIR) -module $(OBJDIR) -c rt.f90 -o $(OBJDIR)/rt.o

$(OBJDIR)/mesh.o: 	$(OBJDIR)/commonRoutines.o | $(OBJDIR)
	$(FC) $(FCFLAGS) -I$(OBJDIR) -module $(OBJDIR) -c mesh.f90 -o $(OBJDIR)/mesh.o

$(OBJDIR)/femHelper.o: 	$(OBJDIR)/commonRoutines.o | $(OBJDIR)
	$(FC) $(FCFLAGS) -I$(OBJDIR) -module $(OBJDIR) -c femHelper.f90 -o $(OBJDIR)/femHelper.o

$(OBJDIR)/rtHelper.o: 	$(OBJDIR)/commonRoutines.o | $(OBJDIR)
	$(FC) $(FCFLAGS) -I$(OBJDIR) -module $(OBJDIR) -c rtHelper.f90 -o $(OBJDIR)/rtHelper.o

$(OBJDIR)/commonRoutines.o: 	| $(OBJDIR)
	$(FC) $(FCFLAGS) -I$(OBJDIR) -module $(OBJDIR) -c commonRoutines.f90 -o $(OBJDIR)/commonRoutines.o

# Clean everything
clean:
	rm -rf $(OBJDIR)
	rm -rf *.out
	rm -rf *.o *.mod

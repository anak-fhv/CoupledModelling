# This makefile uses intel fortran compiler for linux

# Makefile variables
FC      = ifort
# -g: use this flag only if debugging is required
FCFLAGS = -parallel -mkl -traceback -openmp -CB
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
OBJECTS := 	$(OBJDIR)/runModel.o $(OBJDIR)/problem.o \
		   	$(OBJDIR)/mesh.o $(OBJDIR)/fem.o $(OBJDIR)/rt.o \
			$(OBJDIR)/utilities.o

# Make all
all: $(OBJECTS)
	$(FC) $(FCFLAGS) -o $(OBJDIR)/test $(OBJECTS) $(LIBS)

$(OBJDIR)/runModel.o: 	$(OBJDIR)/problem.o \
						$(OBJDIR)/mesh.o $(OBJDIR)/fem.o \
						$(OBJDIR)/rt.o $(OBJDIR)/utilities.o | $(OBJDIR)
	$(FC) $(FCFLAGS) -I$(OBJDIR) -c runModel.f90 -o $(OBJDIR)/runModel.o

$(OBJDIR)/problem.o: 	$(OBJDIR)/rt.o $(OBJDIR)/fem.o \
								$(OBJDIR)/mesh.o $(OBJDIR)/utilities.o | $(OBJDIR)
	$(FC) $(FCFLAGS) -I$(OBJDIR) -module $(OBJDIR) -c problem.f90 -o $(OBJDIR)/problem.o

$(OBJDIR)/mesh.o: 	$(OBJDIR)/utilities.o | $(OBJDIR)
	$(FC) $(FCFLAGS) -I$(OBJDIR) -module $(OBJDIR) -c mesh.f90 -o $(OBJDIR)/mesh.o

$(OBJDIR)/fem.o: 	$(OBJDIR)/utilities.o $(OBJDIR)/mesh.o | $(OBJDIR)
	$(FC) $(FCFLAGS) -I$(OBJDIR) -module $(OBJDIR) -c fem.f90 -o $(OBJDIR)/fem.o

$(OBJDIR)/rt.o: 	$(OBJDIR)/utilities.o $(OBJDIR)/mesh.o | $(OBJDIR)
	$(FC) $(FCFLAGS) -I$(OBJDIR) -module $(OBJDIR) -c rt.f90 -o $(OBJDIR)/rt.o

$(OBJDIR)/utilities.o: 	| $(OBJDIR)
	$(FC) $(FCFLAGS) -I$(OBJDIR) -module $(OBJDIR) -c utilities.f90 -o $(OBJDIR)/utilities.o

# Clean everything
clean:
	rm -rf $(OBJDIR)
	rm -rf *.out
	rm -rf *.o *.mod

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

# Clean everything
clean:
	rm -rf $(OBJDIR)
	rm -rf *.out
	rm -rf *.o *.mod

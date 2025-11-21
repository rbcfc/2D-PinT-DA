#
# Include paths and libraries: change accordingly
# 
#---------------------------------------------------------------------
EXEC_BIN = parareal
NCDF_PATH= -I$(shell brew --prefix netcdf)/include \
		-I$(shell brew --prefix netcdf-fortran)/include
NETCDF_LIB = -L$(shell brew --prefix netcdf)/lib \
		-L$(shell brew --prefix netcdf-fortran)/lib \
		-lnetcdff -lnetcdf
ARPACK_LIB = -L$(shell brew --prefix arpack)/lib -larpack

OPENBLAS_INC = -I$(shell brew --prefix openblas)/include

OPENBLAS_LIB = -L$(shell brew --prefix openblas)/lib -lopenblas
#---------------------------------------------------------------------
F_C = gfortran
F_P =
F_O = -fdefault-real-8 -fdefault-double-8 -fimplicit-none -O3 -g -mtune=native \
 -fopenmp $(NCDF_PATH) $(OPENBLAS_INCLUDE) 
L_O = $(NETCDF_LIB) $(ARPACK_LIB) $(OPENBLAS_LIB) -llapack -lblas
RM  = rm -f
PREF=

OBJSET  = \
	parareal.o shallow_water.o shallow_water_utils.o main.o inexact_cg.o \
	eigenvalues.o shallow_water_step.o shallow_water_adjoint.o shallow_water_step_b.o \
	output.o cpu_time.o dgesvd.o
#
# Target
all: $(EXEC_BIN)
		@echo compilation is OK
# cleaning objects, libraries and executables
clean:
		$(RM) *.o *.mod $(EXEC_BIN)
		@echo \(.o .mod and executables are removed\)


$(EXEC_BIN) :$(OBJSET)
	$(F_C) $(F_O) -o $(EXEC_BIN) $(OBJSET) $(L_O)

%.o: %.f90
		$(F_C) -c $(F_O) $<
%.o: %.F90
		$(F_C) -c $(F_O) $<
%.o: %.f
		$(F_C) -c $(F_O) $<

# dependencies

parareal.o : shallow_water.o shallow_water_adjoint.o eigenvalues.o
shallow_water.o : shallow_water_utils.o shallow_water_step.o eigenvalues.o
shallow_water_adjoint.o : shallow_water.o shallow_water_step_b.o
inexact_cg.o : parareal.o eigenvalues.o
main.o : shallow_water.o inexact_cg.o parareal.o eigenvalues.o output.o shallow_water_adjoint.o 

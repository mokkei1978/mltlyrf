FC 	= ifort
FFLAGS	= -fast
FFLAGS	+= -parallel
#FFLAGS 	= -ipo -O3 -no-prec-div -static
#FFLAGS 	= -axP -vec-report0
#FFLAGS 	= -fast -vec-report0
#FFLAGS 	= -vec-report0
#FFLAGS 	= -check -debug extended -debug-parameters -warn
#FC	= frt
##FFLAGS	= -X9 -fw -Kparallel,instance=2,ocl #-Haesux
#FFLAGS	= -X9 -Am -Haesux
#FC	= gfortran

MODS = 	common.o  \
	euler_asselin.o \
	runge_kutta2.o \
	harmonic_diffusion.o \
	biharmonic_diffusion.o \
	user_module.o

MAIN = 	initial_condition.o \
	diagnostic.o \
	inertial_forcings.o \
	boundary_condition.o \
	monitoring.o \
	output.o \
	available_file_unit.o

.SUFFIXES:

.SUFFIXES: .f .f90 .F90 .o

.f90.o:
	$(FC) $(FFLAGS) $(OPT) -c $< 
.F90.o:
	$(FC) $(FFLAGS) $(OPT) -c $< 
.f.o:
	$(FC) $(FFLAGS) $(OPT) -c $< 


main: $(MODS) $(MAIN) main.o
	$(FC) $(FFLAGS) -o $@.exe $(MODS) $(MAIN) main.o

.PHONY: clean
clean:
	\rm -f *.o *.exe *.mod *~

euler_asselin.o: common.o
runge_kutta2.o: common.o
harmonic_diffusion.o: common.o
biharmonic_diffusion.o: common.o
user_module.o: common.o

initial_condition.o: common.o
diagnostic.o: common.o
inertial_forcings.o: common.o
boundary_condition.o: common.o
monitoring.o: common.o
output.o: common.o

main.o: $(MODS)

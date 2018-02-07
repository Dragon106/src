FC = ifort
FCFLAGS = -Ofast -static -traceback -mkl
#FCFLAGS = -Ofast -static -g -check all -fpe0 -warn -traceback -debug extended
LIBS = 

PROG = HHG_MOD_NEW

all: clean $(PROG) 
$(PROG): hhg.o hhg_mod.o deriv_dip.o splines.o integ_weigs.o read_inputs.o vars.o pars.o utils.o 
	$(FC) $(FCFLAGS) -o $(PROG) hhg.o hhg_mod.o deriv_dip.o splines.o integ_weigs.o read_inputs.o vars.o pars.o utils.o

hhg.o: hhg.f90 hhg_mod.o deriv_dip.o splines.o read_inputs.o vars.o pars.o utils.o
	$(FC) $(FCFLAGS) -c hhg.f90

hhg_mod.o: hhg_mod.f90 integ_weigs.o pars.o vars.o
	$(FC) $(FCFLAGS) -c hhg_mod.f90

deriv_dip.o: deriv_dip.f90
	$(FC) $(FCFLAGS) -c deriv_dip.f90

splines.o: splines.f90
	$(FC) $(FCFLAGS) -c splines.f90

integ_weigs.o: integ_weigs.f90 pars.o
	$(FC) $(FCFLAGS) -c integ_weigs.f90

read_inputs.o: read_inputs.f90 pars.o vars.o
	$(FC) $(FCFLAGS) -c read_inputs.f90

vars.o: vars.f90 pars.o
	$(FC) $(FCFLAGS) -c vars.f90

pars.o: pars.f90
	$(FC) $(FCFLAGS) -c pars.f90

utils.o: utils.f90
	$(FC) $(FCFLAGS) -c utils.f90

clean:
	rm -f *.o *.mod

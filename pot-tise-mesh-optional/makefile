FC      = ifort
FCFLAGS = -Ofast -static -xHOST -parallel
LIBS =
PROG = SAEPOT_TISE_Vc

all: clean $(PROG)

$(PROG): type_vars.o consts.o read_pars.o quadrule.o GAUQUA.o grid.o eledens_info.o edens_tise.o edens_gau.o gradv.o gauss_basis.o SAEPOT.o utils.o
	$(FC) $(FCFLAGS) -o $(PROG) type_vars.o consts.o read_pars.o quadrule.o GAUQUA.o grid.o eledens_info.o edens_tise.o edens_gau.o gradv.o gauss_basis.o SAEPOT.o utils.o

type_vars.o: type_vars.f90
	$(FC) $(FCFLAGS) -c type_vars.f90

consts.o: consts.f90 type_vars.o
	$(FC) $(FCFLAGS) -c consts.f90

read_pars.o: read_pars.f90 type_vars.o
	$(FC) $(FCFLAGS) -c read_pars.f90

quadrule.o: quadrule.f90
	$(FC) $(FCFLAGS) -c quadrule.f90

GAUQUA.o: GAUQUA.f
	$(FC) $(FCFLAGS) -c GAUQUA.f

grid.o: grid.f90 type_vars.o consts.o read_pars.o utils.o
	$(FC) $(FCFLAGS) -c grid.f90

eledens_info.o: eledens_info.f90 type_vars.o consts.o grid.o utils.o
	$(FC) $(FCFLAGS) -c eledens_info.f90

edens_tise.o: edens_tise.f90 type_vars.o consts.o grid.o eledens_info.o utils.o
	$(FC) $(FCFLAGS) -c edens_tise.f90

edens_gau.o: edens_gau.f90 type_vars.o consts.o grid.o eledens_info.o gauss_basis.o utils.o
	$(FC) $(FCFLAGS) -c edens_gau.f90

gradv.o: gradv.f90 grid.o
	$(FC) $(FCFLAGS) -c gradv.f90

gauss_basis.o: gauss_basis.f90 type_vars.o
	$(FC) $(FCFLAGS) -c gauss_basis.f90

SAEPOT.o: SAEPOT.f90 type_vars.o consts.o read_pars.o grid.o eledens_info.o edens_tise.o edens_gau.o gradv.o
	$(FC) $(FCFLAGS) -c SAEPOT.f90

utils.o: utils.f90
	$(FC) $(FCFLAGS) -c utils.f90
	
clean:
	rm -f *.o *.mod

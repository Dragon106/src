FC      = ifort
FCFLAGS = -Ofast -static -xHOST -parallel

PROG = SAEPOT_GAU_Vc

SRC  = consts.f90 type_vars.f90 read_pars.f90 quadrule.f90 GAUQUA.f grid.f90 eledens_info.f90 edens_tise.f90 edens_gau.f90 gradv.f90 gauss_basis.f90 SAEPOT.f90 utils.f90

all: clean $(PROG)

$(PROG): 


clean:
	rm -f *.o *.mod

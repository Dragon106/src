FC = ifort
FCFLAGS = -Ofast -static -xHOST
FCFLAGS += -parallel -openmp 
#FCFLAGS += -opt-report=5 -opt-report-phase=vec,openmp
FCFLAGS += -traceback
SRC = type_vars.f90 consts.f90 quadrule.f90 mesh.f90 laser.f90 PtrArr.f90 data.f90 tdse.f90 coeffs.f90 output.f90 rk4.f90 solve_tdse.f90 main.f90

PROG = TDSE_test

all: $(PROG)
$(PROG): $(SRC) 
	$(FC) $(FCFLAGS) $(SRC) -o $(PROG)
 
clean:
	rm -f *.o *.mod

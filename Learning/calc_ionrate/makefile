FC = ifort
FCFLAGS = -Ofast -static -xHOST -parallel
LIBS = 

SRC  = type_vars.f90 derivfx.f90 calc_rate.f90 main.f90
 
PROG = calc_rate

all: $(PROG) clean

$(PROG): $(SRC)
	$(FC) $(FCFLAGS) $(SRC) -o $(PROG)

clean:
	rm -f *.o *.mod

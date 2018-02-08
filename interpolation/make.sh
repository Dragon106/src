#!/bin/bash

FC='ifort'
FLAGS='-Ofast -static -xHOST -parallel'
LIBS='-mkl -lmkl_lapack95_lp64'
SRC='interpot.f90 cspl.f90 spl1.f90 grid.f90 gradv.f90 read_pars.f90 GAUQUA.f quadrule.f90 consts.f90 type_vars.f90 utils.f90'

$FC -o INTERP $FLAGS $SRC $LIBS
./INTERP

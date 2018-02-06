#!/bin/bash

#source /home/dinhhanh/Tu/Apps/bin/envload

FC='ifort'
FLAGS='-Ofast -static -xHOST -openmp -parallel -mkl'
LIBS='-lmkl_lapack95_lp64'

SRC='basis.f90 GAUQUA.f quadrule.f90 consts.f90 hamil.f90 mesh.f90 molpot.f90 gradv.f90 tise2d.f90 utils.f90 read_pars.f90 wfs.f90' 

$FC $FLAGS -o TISE $SRC $LIBS


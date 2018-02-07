#!/bin/bash

source /home/dinhhanh/Tu/Apps/bin/envload

FC='ifort'
FLAGS='-Ofast -static -xHOST -openmp -mkl=parallel -parallel'

SRC='consts.f90 type_vars.f90 read_pars.f90 quadrule.f90 GAUQUA.f grid.f90 eledens_info.f90 edens_tise.f90 edens_gau.f90 gradv.f90 gauss_basis.f90 SAEPOT.f90 utils.f90' 

$FC $FLAGS -o SAEPOT_GAU_Vc $SRC

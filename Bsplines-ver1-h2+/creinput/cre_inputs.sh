#!/bin/bash 

wdir=RUN

mmin=0
mmax=2
lmax=5

Nrad=400
Neta=100
Rmax=120.0

NBr=200
deg=9
gam=2.0
iline=300

Emax=7.0

mkdir $wdir
cd $wdir

for m in $(seq $mmin $mmax); do 

    mq=$m
#    let lmax=$lmax+$i-$mmin

    mkdir $m
    cd $m
    
    fln='inputs.in'

    echo -e "! MESH PARAMETERS ( NRAD, NETA, RMAX ) \n" > $fln
    
    echo -e "$Nrad  $Neta  $Rmax \n" >> $fln

    echo -e "! BSPLINES PARAMETERS ( NBR, DEG, GAMMA, iline ) \n" >> $fln
        
    echo -e "$NBr  $deg  $gam  $iline \n" >> $fln
    
    echo -e "! ASSOCIATED LEGENDRE POLYNOMIAL PARAMETERS ( LMAX, QM ) \n" >> $fln
    
    echo -e "$lmax $mq \n"  >> $fln

    echo -e "! MAXIMUM ENERGY ( Emax ) \n" >> $fln
    
    echo -e "$Emax \n" >> $fln
    
    echo -e "! POTENTIAL ( POT_FILE_NAME ) \n" >> $fln
    
    echo -e "POT.DAT \n" >> $fln 
      
    cd ..
done 


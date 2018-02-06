#!/bin/bash
#
#$ -N h-b
#$ -l h_rt=500:00:00
#$ -pe orte 24
#$ -cwd
#$ -j y
#$ -S /bin/bash
#

date
echo

wdir=RUN

inp='inputs.in'
pot='POT.DAT'
tise='TISE'

cd $wdir
for d in *; do
             
    echo 
    echo "Running in $d directory"
    echo 
    
    cd $d
    
#    cp ../../$pot . ! for another molecule
    ../../$tise $inp > out.log
 
    cd ..

    echo     
    echo "Finish in $d and continue..."
    echo
 
done 
cd ..

echo
date

#!/bin/bash 

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
    
#    cp ../../$pot .
    ../../$tise $inp > out.log
 
    cd ..

    echo     
    echo "Finish in $d and continue..."
    echo
 
done 
cd ..


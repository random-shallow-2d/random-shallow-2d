#!/bin/bash

#L = 1D size
L_init=$1
k=$2
#m = number of iterations
m=$3
#numCircuits = number of circuits to do
numCircuits=$4
#numReps = how many times to run a circuit
numReps=$5



#method= 1 for bond and 2 for err
method=$6
#err=truncation error
err=$7
#maxB
maxB=$8
directory="results/brickwork/L${L_init}_k${k}_m${m}_numCircuits${numCircuits}_numReps${numReps}_method${method}_err${err}_maxB${maxB}"
mkdir -p "$directory"
for (( q=0; q<=$k-1; q++))
    do
    L=$(($L_init+8*$q))
    num_unit=$(( (( (($L-1)) + 2 * (($L-1)) / 8)) * $m ))
    for (( i=1; i<=$numCircuits; i++ ))
        do
            #filename="/Volumes/external/circuits/brickwork_L9-401_550its_50reps/circuit$i"
            filename="$directory/circuit$i"
            filename+="_L$L"
            filename+="_m$m"
            filename+=".csv"
            python matgen.py $num_unit $filename

            for (( j=1; j<=$numReps; j++ ))
                do
                    outfile="$directory/result_circuit$i"
                    outfile+="_L$L"
                    outfile+="_m$m"
                    outfile+="_rep$j"
                    outfile+=".txt"
                    ./brickwork $method $err $maxB $L $filename $m $outfile
                done
            rm "$filename"
        done
    done

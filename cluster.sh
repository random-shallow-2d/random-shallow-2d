#!/bin/bash

#L = 1D size
L_low=$1
L_high=$2
#m = number of iterations
m=$3
#numCircuits = number of circuits to do
numCircuits=$4


#method= 1 for bond and 2 for err
method=$5
#err=truncation error
err=$6
#maxB
maxB=$7
increment=$8
directory="results/cluster/Llow${L_low}_Lhigh${L_high}_m${m}_numCircuits${numCircuits}_method${method}_err${err}_maxB${maxB}_increment${increment}"
mkdir -p "$directory"
for (( L=$L_low; L<=$L_high; L=L+$increment ))
do
	for (( i=1; i<=$numCircuits; i++ ))
    do
        outfile="${directory}/result_circuit$i"
		outfile+="_L$L"
        outfile+="_m$m"
        outfile+=".txt"
        ./cluster $method $err $maxB $L $m $outfile
    done
done

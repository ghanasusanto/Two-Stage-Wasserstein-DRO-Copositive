#!/bin/bash
cd /vol/bitbucket/gah11/two_stage/vs_decision_rules

#/usr/bin/matlab -nosplash -nojvm -nodesktop -nodisplay -r "run"


for (( d_=6; d_<=7; d_++ ))
do
for (( N_=1; N_<=8; N_++ ))
do
for (( i=1; i<=100; i++ ))
do
FILE="results/${d_}_${N_}/d_${d_}_N_${N_}_i_${i}.mat"
FILE2="runs/d_${d_}_N_${N_}_i_${i}"
if  [ ! -f $FILE ] && [ ! -f $FILE2 ]
then
    touch $FILE2
    echo $d_ $N_ $i
    /usr/bin/matlab -nosplash -nojvm -nodesktop -nodisplay -r "run_d_N_quit(${d_},${N_},${i})"
fi
done
done
done

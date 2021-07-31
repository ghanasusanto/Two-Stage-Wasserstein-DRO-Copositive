#!/bin/bash
echo $1
for (( N_=1; N_<=8; N_++ ))
do
for (( i=1; i<=100; i++ ))
do
FILE="results_new/${1}_${N_}/d_${1}_N_${N_}_i_${i}.mat"
if  [ ! -f $FILE ]
then
    echo $1 $N_ $i
    echo -e "universe = vanilla \n Notification = Error \n environment = MOSEKLM_LICENSE_FILE=/homes/gah11/mosek/mosek.lic  \n executable      = /usr/bin/matlab \n Arguments	= \"-nosplash -nojvm -nodesktop -nodisplay -r 'run_d_N_new($1,$N_,$i)'\" \n \n Requirements= Arch == \"x86_64\" && TotalMemory > 12000 \n queue 1 " > commands/trial_${1}_${N_}_${i}.cmd
    condor_submit commands/trial_${1}_${N_}_${i}.cmd
fi
done
done


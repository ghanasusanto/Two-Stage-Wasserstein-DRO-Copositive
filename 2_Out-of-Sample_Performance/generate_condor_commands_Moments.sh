#!/bin/bash
rm
echo $1
for (( i=1; i<=100; i++ ))
do
FILE="results/${1}/N_${1}_i_${i}_Moments.mat"
if  [ ! -f $FILE ]
then
    echo $1 $i
    echo -e "universe = vanilla \n Notification = Error \n environment = MOSEKLM_LICENSE_FILE=/homes/gah11/mosek/mosek.lic  \n executable      = /usr/bin/matlab \n Arguments	= \"-nosplash -nojvm -nodesktop -nodisplay -r 'run_moments_N_i($1,$i)'\" \n \n Requirements= Arch == \"x86_64\" && TotalMemory > 12000 \n queue 1 " > commands/trial_${1}_${i}.cmd
    condor_submit commands/trial_${1}_${i}.cmd
fi
done

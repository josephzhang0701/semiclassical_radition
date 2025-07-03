#!/bin/bash

# this script is used to subsitute the parameters in "qsub.txt", and then generate the qsub script file

loop=1
step=1
echo
echo "qsub_script_generation started"
for particle_nr in `seq 1 $step $loop`
do
  sed -e "s/particle_nr/$particle_nr/g" ./qjob2.txt >./qjobs/qjob$particle_nr.sh
  chmod u+x ./qjobs/qjob$particle_nr.sh
done
echo "qsub_script_generation finished"
echo
exit 0 

#!/bin/bash
./qjobs_generator.sh
loop=1
step=1
echo
echo "submitting jobs"
for particle_nr in `seq 1 $step $loop`
do
  qsub ./qjobs/qjob$particle_nr.sh
done
echo "all jobs have been submitted"
echo

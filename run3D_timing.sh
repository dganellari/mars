#!/bin/bash

SP_PATH=/Users/liudmilakaragyaur/Documents/Eurohack2020_local/mars/build/st_example

for level in {5,10,20,30,40,50}
do
   srun -n4 --ntasks-per-core=1 --ntasks-per-node=1 --cpu-bind=no $SP_PATH $level 0
done 
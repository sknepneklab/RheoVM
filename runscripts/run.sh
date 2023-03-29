#!/bin/bash

num_period=(40 40 40 20 20 20 20 20 10 10 10 10 10 6 6 6)
for((i=0;i<$((15));i++))
do
    python run_shear_frequency_sweep.py ${i} ${num_period[i]}
done

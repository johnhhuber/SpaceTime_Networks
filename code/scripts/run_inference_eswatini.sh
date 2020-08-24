#!/bin/bash
#$ -pe smp 8
#$ -q long
#$ -N eswatini_inference
#$ -t 1-5

module load gcc
module load gsl

# run inferences
inDir=../../data
outDir=../../data/eswatini_inference

../src/main.out $(inDir)/eswatini_inference/nodes.csv $(inDir)/network_sim/sources_NULL.csv $(inDir)/network_sim/settings_space_time_travel.csv $(ouDir)/space_time_travel/replicate_$SGE_TASK_ID/networks.csv $(outDir)/space_time_travel/replicate_$SGE_TASK_ID/scalars.csv $(outDir)/space_time_travel/replicate_$SGE_TASK_ID/swaps.csv $SGE_TASK_ID
../src/main.out $(inDir)/eswatini_inference/nodes.csv $(inDir)/network_sim/sources_NULL.csv $(inDir)/network_sim/settings_space_time_bth.csv $(ouDir)/space_time_bth/replicate_$SGE_TASK_ID/networks.csv $(outDir)/space_time_bth/replicate_$SGE_TASK_ID/scalars.csv $(outDir)/space_time_bth/replicate_$SGE_TASK_ID/swaps.csv $SGE_TASK_ID
../src/main.out $(inDir)/eswatini_inference/nodes.csv $(inDir)/network_sim/sources_NULL.csv $(inDir)/network_sim/settings_space_time.csv $(ouDir)/space_time/replicate_$SGE_TASK_ID/networks.csv $(outDir)/space_time/replicate_$SGE_TASK_ID/scalars.csv $(outDir)/space_time/replicate_$SGE_TASK_ID/swaps.csv $SGE_TASK_ID
../src/main.out $(inDir)/eswatini_inference/nodes.csv $(inDir)/network_sim/sources_NULL.csv $(inDir)/network_sim/settings_time_travel.csv $(ouDir)/time_travel/replicate_$SGE_TASK_ID/networks.csv $(outDir)/time_travel/replicate_$SGE_TASK_ID/scalars.csv $(outDir)/time_travel/replicate_$SGE_TASK_ID/swaps.csv $SGE_TASK_ID
../src/main.out $(inDir)/eswatini_inference/nodes.csv $(inDir)/network_sim/sources_NULL.csv $(inDir)/network_sim/settings_time_bth.csv $(ouDir)/time_bth/replicate_$SGE_TASK_ID/networks.csv $(outDir)/time_bth/replicate_$SGE_TASK_ID/scalars.csv $(outDir)/time_bth/replicate_$SGE_TASK_ID/swaps.csv $SGE_TASK_ID

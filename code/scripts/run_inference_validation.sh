#!/bin/bash
#$ -pe smp 8
#$ -q long
#$ -N eswatini_inference
#$ -t 1-5

module load gcc
module load gsl

# run inferences
inDir=../../data/
outDir=../../data/eswatini_inference

../src/main.out $(outDir)/space_time_travel/validation/nodes_0_FULL.csv $(inDir)/network_sim/sources_NULL.csv $(inDir)/network_sim/settings_space_time_travel.csv $(ouDir)/space_time_travel/validation/replicate_$SGE_TASK_ID/networks_0.csv $(outDir)/space_time_travel/validation/replicate_$SGE_TASK_ID/scalars_0.csv $(outDir)/space_time_travel/validation/replicate_$SGE_TASK_ID/swaps_0.csv $SGE_TASK_ID
../src/main.out $(outDir)/space_time_bth/validation/nodes_0_FULL.csv $(inDir)/network_sim/sources_NULL.csv $(inDir)/network_sim/settings_space_time_bth.csv $(ouDir)/space_time_bth/validation/replicate_$SGE_TASK_ID/networks_0.csv $(outDir)/space_time_bth/validation/replicate_$SGE_TASK_ID/scalars_0.csv $(outDir)/space_time_bth/validation/replicate_$SGE_TASK_ID/swaps_0.csv $SGE_TASK_ID
../src/main.out $(outDir)/space_time/validation/nodes_0_FULL.csv $(inDir)/network_sim/sources_NULL.csv $(inDir)/network_sim/settings_space_time.csv $(ouDir)/space_time/validation/replicate_$SGE_TASK_ID/networks_0.csv $(outDir)/space_time/validation/replicate_$SGE_TASK_ID/scalars_0.csv $(outDir)/space_time/validation/replicate_$SGE_TASK_ID/swaps_0.csv $SGE_TASK_ID
../src/main.out $(outDir)/time_travel/validation/nodes_0_FULL.csv $(inDir)/network_sim/sources_NULL.csv $(inDir)/network_sim/settings_time_travel.csv $(ouDir)/time_travel/validation/replicate_$SGE_TASK_ID/networks_0.csv $(outDir)/time_travel/validation/replicate_$SGE_TASK_ID/scalars_0.csv $(outDir)/time_travel/validation/replicate_$SGE_TASK_ID/swaps_0.csv $SGE_TASK_ID
../src/main.out $(outDir)/time_bth/validation/nodes_0_FULL.csv $(inDir)/network_sim/sources_NULL.csv $(inDir)/network_sim/settings_time_bth.csv $(ouDir)/time_bth/validation/replicate_$SGE_TASK_ID/networks_0.csv $(outDir)/time_bth/validation/replicate_$SGE_TASK_ID/scalars_0.csv $(outDir)/time_bth/validation/replicate_$SGE_TASK_ID/swaps_0.csv $SGE_TASK_ID

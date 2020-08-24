#!/bin/bash
#$ -q long
#$ -N sweep
#$ -t 1-2000

module load R

# generate data
Rscript run_gen_data_sweep.R $SGE_TASK_ID

# run inferences
inDir=../../data/sim_sweep/data
outDir=../../output/sim_sweep

../src/main.out $(inDir)/nodes_$((SGE_TASK_ID-1))_FULL.csv ../../data/network_sim/sources_NULL.csv ../../data/network_sim/settings_space_time_travel.csv $(outDir)/networks/networks_stt_$((SGE_TASK_ID-1)).csv $(outDir)/scalars/scalars_stt_$((SGE_TASK_ID-1)).csv $(outDir)/swaps/swaps_stt_$((SGE_TASK_ID-1)).csv $SGE_TASK_ID
../src/main.out $(inDir)/nodes_$((SGE_TASK_ID-1))_FULL.csv ../../data/network_sim/sources_NULL.csv ../../data/network_sim/settings_space_time_bth.csv $(outDir)/networks/networks_stb_$((SGE_TASK_ID-1)).csv $(outDir)/scalars/scalars_stb_$((SGE_TASK_ID-1)).csv $(outDir)/swaps/swaps_stb_$((SGE_TASK_ID-1)).csv $SGE_TASK_ID
../src/main.out $(inDir)/nodes_$((SGE_TASK_ID-1))_FULL.csv ../../data/network_sim/sources_NULL.csv ../../data/network_sim/settings_space_time.csv $(outDir)/networks/networks_stn_$((SGE_TASK_ID-1)).csv $(outDir)/scalars/scalars_stn_$((SGE_TASK_ID-1)).csv $(outDir)/swaps/swaps_stn_$((SGE_TASK_ID-1)).csv $SGE_TASK_ID

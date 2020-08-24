#!/bin/bash
#$ -q long
#$ -N process_sweep
#$ -t 1-2000

module load R

# compute accuracy metrics
Rscript process_accuracy_sweep.R ../../output/sim_sweep/networks/ ../../data/sim_sweep/data/ ../../output/sim_sweep/accuracy/ $SGE_TASK_ID stt
Rscript process_accuracy_sweep.R ../../output/sim_sweep/networks/ ../../data/sim_sweep/data/ ../../output/sim_sweep/accuracy/ $SGE_TASK_ID stb
Rscript process_accuracy_sweep.R ../../output/sim_sweep/networks/ ../../data/sim_sweep/data/ ../../output/sim_sweep/accuracy/ $SGE_TASK_ID stn

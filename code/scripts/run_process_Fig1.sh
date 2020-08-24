#!/bin/bash
#$ -N process_fig_1
#$ -q long

module load R

# process MC3 for convergence
Rscript ../R/process_MC3.R ../../data/figs/fig_1/10local_90imported/space_time_travel/
Rscript ../R/process_MC3.R ../../data/figs/fig_1/10local_90imported/space_time_bth/
Rscript ../R/process_MC3.R ../../data/figs/fig_1/10local_90imported/space_time/
Rscript ../R/process_MC3.R ../../data/figs/fig_1/10local_90imported/time_travel/
Rscript ../R/process_MC3.R ../../data/figs/fig_1/10local_90imported/time_bth/

Rscript ../R/process_MC3.R ../../data/figs/fig_1/50local_50imported/space_time_travel/
Rscript ../R/process_MC3.R ../../data/figs/fig_1/50local_50imported/space_time_bth/
Rscript ../R/process_MC3.R ../../data/figs/fig_1/50local_50imported/space_time/
Rscript ../R/process_MC3.R ../../data/figs/fig_1/50local_50imported/time_travel/
Rscript ../R/process_MC3.R ../../data/figs/fig_1/50local_50imported/time_bth/

Rscript ../R/process_MC3.R ../../data/figs/fig_1/95local_05imported/space_time_travel/
Rscript ../R/process_MC3.R ../../data/figs/fig_1/95local_05imported/space_time_bth/
Rscript ../R/process_MC3.R ../../data/figs/fig_1/95local_05imported/space_time/
Rscript ../R/process_MC3.R ../../data/figs/fig_1/95local_05imported/time_travel/
Rscript ../R/process_MC3.R ../../data/figs/fig_1/95local_05imported/time_bth/

# compute accuracy metrics
Rscript ../R/process_accuracy.R ../../data/figs/fig_1/10local_90imported/space_time_travel ../../data/figs/fig_1/10local_90imported/network_true_1_FULL.csv
Rscript ../R/process_accuracy.R ../../data/figs/fig_1/10local_90imported/space_time_bth ../../data/figs/fig_1/10local_90imported/network_true_1_FULL.csv
Rscript ../R/process_accuracy.R ../../data/figs/fig_1/10local_90imported/space_time ../../data/figs/fig_1/10local_90imported/network_true_1_FULL.csv
Rscript ../R/process_accuracy.R ../../data/figs/fig_1/10local_90imported/time_travel ../../data/figs/fig_1/10local_90imported/network_true_1_FULL.csv
Rscript ../R/process_accuracy.R ../../data/figs/fig_1/10local_90imported/time_bth ../../data/figs/fig_1/10local_90imported/network_true_1_FULL.csv

Rscript ../R/process_accuracy.R ../../data/figs/fig_1/50local_50imported/space_time_travel ../../data/figs/fig_1/50local_50imported/network_true_1_FULL.csv
Rscript ../R/process_accuracy.R ../../data/figs/fig_1/50local_50imported/space_time_bth ../../data/figs/fig_1/50local_50imported/network_true_1_FULL.csv
Rscript ../R/process_accuracy.R ../../data/figs/fig_1/50local_50imported/space_time ../../data/figs/fig_1/50local_50imported/network_true_1_FULL.csv
Rscript ../R/process_accuracy.R ../../data/figs/fig_1/50local_50imported/time_travel ../../data/figs/fig_1/50local_50imported/network_true_1_FULL.csv
Rscript ../R/process_accuracy.R ../../data/figs/fig_1/50local_50imported/time_bth ../../data/figs/fig_1/50local_50imported/network_true_1_FULL.csv

Rscript ../R/process_accuracy.R ../../data/figs/fig_1/95local_05imported/space_time_travel ../../data/figs/fig_1/95local_05imported/network_true_1_FULL.csv
Rscript ../R/process_accuracy.R ../../data/figs/fig_1/95local_05imported/space_time_bth ../../data/figs/fig_1/95local_05imported/network_true_1_FULL.csv
Rscript ../R/process_accuracy.R ../../data/figs/fig_1/95local_05imported/space_time ../../data/figs/fig_1/95local_05imported/network_true_1_FULL.csv
Rscript ../R/process_accuracy.R ../../data/figs/fig_1/95local_05imported/time_travel ../../data/figs/fig_1/95local_05imported/network_true_1_FULL.csv
Rscript ../R/process_accuracy.R ../../data/figs/fig_1/95local_05imported/time_bth ../../data/figs/fig_1/95local_05imported/network_true_1_FULL.csv

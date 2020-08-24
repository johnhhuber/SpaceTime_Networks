#!/bin/bash
#$ -q long
#$ -N process_validation

module load R

# check for convergence
Rscript process_covergence.R ../../data/eswatini_inference/space_time_travel/validation
Rscript process_covergence.R ../../data/eswatini_inference/space_time_bth/validation
Rscript process_covergence.R ../../data/eswatini_inference/space_time/validation
Rscript process_covergence.R ../../data/eswatini_inference/time_travel/validation
Rscript process_covergence.R ../../data/eswatini_inference/time_bth/validation

# computate accuracy metrics
Rscript process_accuracy.R ../../data/eswatini_inference/space_time_travel/validation ../../data/eswatini_inference/space_time_travel/validation/network_true_0_FULL.csv
Rscript process_accuracy.R ../../data/eswatini_inference/space_time_bth/validation ../../data/eswatini_inference/space_time_bth/validation/network_true_0_FULL.csv
Rscript process_accuracy.R ../../data/eswatini_inference/space_time/validation ../../data/eswatini_inference/space_time/validation/network_true_0_FULL.csv
Rscript process_accuracy.R ../../data/eswatini_inference/time_travel/validation ../../data/eswatini_inference/time_travel/validation/network_true_0_FULL.csv
Rscript process_accuracy.R ../../data/eswatini_inference/time_bth/validation ../../data/eswatini_inference/time_bth/validation/network_true_0_FULL.csv

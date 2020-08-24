#!/bin/bash
#$ -q long 
#$ -N process_eswatini

module load R

Rscript process_covergence.R ../../data/eswatini_inference/space_time_travel
Rscript process_covergence.R ../../data/eswatini_inference/space_time_bth
Rscript process_covergence.R ../../data/eswatini_inference/space_time
Rscript process_covergence.R ../../data/eswatini_inference/time_travel
Rscript process_covergence.R ../../data/eswatini_inference/time_bth

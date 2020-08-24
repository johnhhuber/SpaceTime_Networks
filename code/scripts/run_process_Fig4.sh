#!/bin/bash
#$ -q long 
#$ -N process_fig4

module load R 

Rscript ../R/run_Fig4.R

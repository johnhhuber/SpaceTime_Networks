# Inferring person-to-person networks of pathogen transmission: is routine surveillance up to the task?
## John H. Huber, Michelle S. Hsiang, Nomcebo Dlamini, Maxwell Murphy, Sibonakaliso Vilakati, Nomcebo Nhlabathi, Anita Lerch, Rasmus Nielsen, Nyasatu Ntshalintshali, Bryan Greenhouse, T. Alex Perkins

A brief description of the code and data used to run the analyses of Huber et al. (2020) are provided below. 

## Getting started

The scripts for the R and C++ code are written with relative paths, assuming you have the following folder structure:
```
SpaceTimeNetworks_Master
|   README.md
└─── code 
     └─── R
     └─── scripts
     └─── src	 	
└─── data
└─── output    	 
```

All C++ code for the transmission network inference algorithm can be found in code/src/. Bash scripts to run the network inference algorithm on various data sets and process the corresponding outputs can be found in code/scripts/. All R scripts to process inference results and generate figures can be found in code/R/. 

### Software and Package Requirements

The transmission network inference algorithm requires:

* GCC: GNU Compiler Collection (We used gcc/8.3.0)
* GSL: GNU Scientific Library (We used gsl/gcc/2.5)

Various makefiles are provided in code/src/ to build the software project. You may need to adjust the 'g++' command to match your system requirement. 

The R scripts require the following packages:

* boot
* coda 
* grDevices
* igraph
* maptools
* MASS
* mgcv
* pomp
* raster
* rcartcolor
* RColorBrewer
* seqinr
* spatstat 
* tools 

The R scripts should automatically install the necessary packages if they are not already installed. 

## Analysis 

### 1. Run the network inference algorithm on the ideal cases. 

* The script to run the network inference algorithm is available at 'code/scripts/run_inference_Fig1.sh'. 
* After that is complete, the results can be processed and the accuracy metrics can be computed using 'code/scripts/run_process_Fig1.sh'. 
* Finally, the figure can be generated with 'code/R/fig_Fig1.R'.

### 2. Run the network inference algorithm on the Eswatini data.

* The script to run the network inference algorithm on the Eswatini data set is available at 'code/scripts/run_inference_eswatini.sh'. 
* After that is complete, the results can be processed using 'code/scripts/run_process_eswatini.sh'. 
* The figures can be generated using the files 'code/R/fig_Fig2.R' and 'code/R/fig_Fig3.R'. 
* For Figure 4, please run 'code/scripts/run_process_Fig4.sh' and then corresponding R script 'code/R/fig_Fig4.R'. 

### 3. Validate the inferences on the Eswatini data 

* Code to generate the validation data sets can be found at 'code/R/run_validate_eswatini_inferences.R'. These files, however, are already generated. 
* To run the inference algorithm on the existing files, please use 'code/scripts/run_inference_validation.sh'. 
* The results can then be processed and the accuracy metrics can be computed using 'code/scripts/run_process_validation.sh'. 
* Finally, the figures can be generated using 'code/R/fig_Fig5.R' and 'code/R/fig_Fig6.R'.  

### 4. Run the simulation sweep

* Code to run the simulation sweep can be found at 'code/scripts/run_inference_sweep.sh'. 
* Results can be processed using 'code/scripts/run_process_sweep.sh'. 
* To generate figures, please run 'code/R/run_Fig7.R' and 'code/R/fig_Fig7.R' as well as 'code/R/fig_FigS8.R', 'code/R/fig_FigS9.R', and 'code/R/fig_FigS10.R'. 

## Contact Me 

If you encounter difficulties running the code or need help adapting it to meet your own needs, please do not hesitate to contact me. I can be reached at jhuber3 AT nd DOT edu 




#!/bin/bash
#SBATCH -J MCC
#SBATCH -N 1
#SBATCH -n 28
#SBATCH -o %N-%j.log
#SBATCH -e %N-%j.err
#SBATCH -p compute

module add R/3.4.1 gcc/6.3.0

R CMD BATCH Unstratified_BSM_BEAST_MCC_Viper.R

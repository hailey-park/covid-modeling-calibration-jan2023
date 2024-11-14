###################################################################################################
#Title: Run mcmc with different initial seroprevalence
#Author: Hailey Park
#Date: July 25, 2024
###################################################################################################


rm(list=ls())
gc()

#Load libraries
library(dplyr)
library(ggplot2)
library(tibble)
library(lubridate)
library(here)
library(data.table)
library(tidyverse)
library(deSolve) 
library(reshape2)
library(dtw)
library(zoo)

options(dplyr.summarise.inform = FALSE)

# Observed data
observed_data <- read.csv("data/clean-data/severe_cases_inc_extended_oct2024.csv", check.names = FALSE)[,-1] 

# Reading model scripts
source(here::here("model-pop-init-12R-newStart-jan2023.R"))

#Simulation vars
time <- seq(0, 657, by=1) # Number of days. (January 1, 2023 - October 19, 2024)
latent_period <- 3        # Average latent period length.
infectious_period <- 5    # Average infectious period length.
gamma <- 1/infectious_period       
delta <- 1/latent_period
N <- 39512223             # Total population size


#Set simulation size
sim_size <- 10000

#Set up folder structure to save simulation results
dir.create("simulation-results")

args <- commandArgs(trailingOnly = TRUE)
job_selector <- (as.numeric(args[1])-1)

###################################################################################################

if (job_selector==0){ 
  
  #Read time-since files (Jan 1, 2022)
  combined_time_since_by_age <- read.csv("data/clean-data/combined_time_since_distributions_by_age_newDate_jan2023.csv")[,-1]
  
  #Read mcmc functions
  source(here::here("mcmc-calibration-functions-dtw-novelVar-strict-extended-jan2023.R"))
  source(here::here("seir-model-12R-novelVar-extended-jan2023.R"))

  #Set initial params
  starting_params <- c(severity_0_17 = 0.0002487624,
                       severity_18_49 = 0.0020000000,
                       severity_50_64 = 0.0063188657,
                       severity_65_74 = 0.0216498253,
                       severity_75plus = 0.1000000000,
                       nonsevere_waning_rate_vaccine = 5.97682,
                       nonsevere_waning_rate_hybrid = 2.87,
                       lambda = 0.71,
                       severe_waning_level_vaccine = 0.028,
                       severe_waning_level_hybrid = 0.0907,
                       lambda_BA.1.1 = 0.6,
                       lambda_BA.2 = 0.6,
                       lambda_BA.4_5 = 0.6,
                       lambda_BQ = 0.6,
                       lambda_XBB = 0.6,
                       lambda_EG_HV = 0.6,
                       lambda_JN = 0.8,
                       lambda_KP_LB = 0.8)
  #Run MCMC
  chain = run_metropolis_MCMC(starting_params, sim_size)
  
  #Save results
  write.csv(chain, "simulation-results/mcmc-chain-10000sims-test25-novelVar-strict-extended-jan2023.csv")
  
  
} else if (job_selector==1){ #Updated acceptance proposal function, less strict
  
  #Read time-since files (Jan 1, 2022)
  combined_time_since_by_age <- read.csv("data/clean-data/combined_time_since_distributions_by_age_newDate_jan2023.csv")[,-1]
  
  #Read mcmc functions
  source(here::here("mcmc-calibration-functions-dtw-novelVar-lessStrict-extended-jan2023.R"))
  source(here::here("seir-model-12R-novelVar-extended-jan2023.R"))
  
  #Set initial params
  starting_params <- c(severity_0_17 = 0.0002487624,
                       severity_18_49 = 0.0020000000,
                       severity_50_64 = 0.0063188657,
                       severity_65_74 = 0.0216498253,
                       severity_75plus = 0.1000000000,
                       nonsevere_waning_rate_vaccine = 5.97682,
                       nonsevere_waning_rate_hybrid = 2.87,
                       lambda = 0.71,
                       severe_waning_level_vaccine = 0.028,
                       severe_waning_level_hybrid = 0.0907,
                       lambda_BA.1.1 = 0.6,
                       lambda_BA.2 = 0.6,
                       lambda_BA.4_5 = 0.6,
                       lambda_BQ = 0.6,
                       lambda_XBB = 0.6,
                       lambda_EG_HV = 0.6,
                       lambda_JN = 0.8,
                       lambda_KP_LB = 0.8)
  
  #Run MCMC
  chain = run_metropolis_MCMC(starting_params, sim_size)
  
  #Save results
  write.csv(chain, "simulation-results/mcmc-chain-10000sims-test26-novelVar-lessStrict-extended-jan2023.csv")
  
  
} 




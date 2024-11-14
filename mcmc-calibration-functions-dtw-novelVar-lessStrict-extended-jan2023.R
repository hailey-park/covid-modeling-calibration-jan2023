########################################################################################################################
#Title: MCMC Calibration
#Author: Hailey Park
#Date: April 16th, 2024
########################################################################################################################

prediction <- function (initial_conditions, params, contact_matrix_adj) {

  #Incorporate the contact matrix factors into the set of params
  add_params <- c(params, contact_matrix_adj)

  #Run model
  results <- lsoda(initial_conditions, time, seir_model, add_params)
  
  #Reformat results into weekly incidence estimates by age group (separate columns)
  reformatted_results <- as.data.frame(results) %>% dplyr::select(time, new_severe_cases1:new_severe_cases5) %>%
    mutate(week = time %/% 7) %>% group_by(week) %>% summarise(across(c(new_severe_cases1:new_severe_cases5), sum)) %>%
    mutate(age_0_17_inc = new_severe_cases1/total_pop_by_age$total_pop[1] * 100000,
           age_18_49_inc = new_severe_cases2/total_pop_by_age$total_pop[2] * 100000,
           age_50_64_inc = new_severe_cases3/total_pop_by_age$total_pop[3] * 100000,
           age_65_74_inc = new_severe_cases4/total_pop_by_age$total_pop[4] * 100000,
           age_75plus_inc = new_severe_cases5/total_pop_by_age$total_pop[5] * 100000) %>%
    dplyr::select(-c(new_severe_cases1:new_severe_cases5))

  return(reformatted_results)
  #return(results)
}


# The score is the DTW distance of the simulated inc predictions from the observed data.
score <- function (sim_pred, data) {
  data <- data[53:147, ]
  age_0_17 <- dtw(data$`0-17 years`, sim_pred$age_0_17_inc,
                  window.type = "sakoechiba",
                  window.size = 2,
                  keep=TRUE)$distance
  age_18_49 <- dtw(data$`18-49 years`, sim_pred$age_18_49_inc,
                   window.type = "sakoechiba",
                   window.size = 2,
                   keep=TRUE)$distance
  age_50_64 <- dtw(data$`50-64 years`, sim_pred$age_50_64_inc,
                   window.type = "sakoechiba",
                   window.size = 2,
                   keep=TRUE)$distance
  age_65_74 <- dtw(data$`65-74 years`, sim_pred$age_65_74_inc,
                   window.type = "sakoechiba",
                   window.size = 2,
                   keep=TRUE)$distance
  age_75_plus <- dtw(data$`75+ years`, sim_pred$age_75plus_inc,
                     window.type = "sakoechiba",
                     window.size = 2,
                     keep=TRUE)$distance
  
  total <- sum(age_0_17, age_18_49, age_50_64, age_65_74, age_75_plus)
  return(total)
}

sigmoid <- function(x, midpoint = -20, max_y = 10000, min_y = 0, max_diff_x = 200, max_diff_y = 300){
  k <- ((max_y - max_diff_y)/(max_diff_y - min_y))^(1/(max_diff_x - midpoint))
  y <- (max_y * k^midpoint + min_y * k^x) / (k^midpoint  + k^x)
  y_rescaled <- y/10000
  return(y_rescaled)
}

# The likelihood is the SSE scores.
likelihood <- function (params) {
  
  #Parameters
  severity_by_age <- c(params[1], params[2], params[3], params[4], params[5])
  nonsevere_waning_rate_vaccine <- params[6]
  nonsevere_waning_rate_hybrid <- params[7]
  lambda <- params[8]
  severe_waning_level_vacc <- params[9]
  severe_waning_level_hybrid <- params[10]

  #Model initialization
  initial_conditions <- model_pop_init(severity_by_age)
  
  #If initial conditions have negative totals in the susceptible compartments (if severity parameters are too small),
  # return null
  if(any(is.na(initial_conditions))) {
    return(Inf)
  }
  
  #If baseline case-hospitalization fractions by age are not monotonically increasing, 
  # return null
  if(!all(severity_by_age == cummax(severity_by_age))){
    return(Inf)
  }
  
  #Calculate contact matrix adjustment factors
  contact_matrix_adj <- contact_matrix_init(initial_conditions)
  
  #Get predictions
  sim_pred <- prediction(initial_conditions, params, contact_matrix_adj)
  
  
  return((score(sim_pred, observed_data)))
  
}


posterior <- function(param){
  likel <- likelihood(param)
  #priors <- prior(param)
  #print(paste0("Likelihood: ", likel))
  #print(paste0("Priors: ", priors))
  return (likel)# + priors)
}


# Choosing a new parameter value by sampling from a multivariate normal distribution
# centered at the current value, that is called the proposal function. The SDs 
# roughly correspond to the step-size of the chain for each parameter.
proposalfunction <- function(param){
  
  severity_0_17 <- min(max(rnorm(1, mean=param[1], sd=0.000005), 0.000005), 0.0005)
  severity_18_49 <- min(max(rnorm(1, mean=param[2], sd=0.00005),0.0001), 0.002)
  severity_50_64 <- min(max(rnorm(1, mean=param[3], sd=0.0005), 0.001), 0.01)
  severity_65_74 <- min(max(rnorm(1, mean=param[4], sd=0.005), 0.005), 0.1)
  severity_75plus <- min(max(rnorm(1, mean=param[5], sd=0.005), 0.005), 0.1)
  pe_nonsevere_rate_vacc = min(max(rnorm(1, mean=param[6], sd=0.1), 0.3), 6)
  pe_nonsevere_rate_hybrid = min(max(rnorm(1, mean=param[7], sd=0.1), 0.3), 6)
  lambda = min(max(rnorm(1, mean=param[8], sd=0.01), 0.3), 2)
  pe_severe_level_vacc = min(max(rnorm(1, mean=param[9], sd=0.01), 0), 0.1)
  pe_severe_level_hybrid = min(max(rnorm(1, mean=param[10], sd=0.01), 0), 0.1)
  lambda_BA.1.1 = min(max(rnorm(1, mean=param[11], sd=0.01), 0.3), 1.2)
  lambda_BA.2 = min(max(rnorm(1, mean=param[12], sd=0.01), 0.3), 1.2)
  lambda_BA.4_5 = min(max(rnorm(1, mean=param[13], sd=0.01), 0.3), 1.2)
  lambda_BQ = min(max(rnorm(1, mean=param[14], sd=0.01), 0.3), 1.2)
  lambda_XBB = min(max(rnorm(1, mean=param[15], sd=0.01), 0.3), 1.2)
  lambda_EG_HV = min(max(rnorm(1, mean=param[16], sd=0.01), 0.3), 1.2)
  lambda_JN = min(max(rnorm(1, mean=param[17], sd=0.01), 0.3), 1.2)
  lambda_KP_LB = min(max(rnorm(1, mean=param[18], sd=0.01), 0.3), 1.2)
  
  return(c(severity_0_17, severity_18_49, severity_50_64, severity_65_74, severity_75plus,
           pe_nonsevere_rate_vacc, pe_nonsevere_rate_hybrid, lambda, pe_severe_level_vacc, pe_severe_level_hybrid,
           lambda_BA.1.1, lambda_BA.2, lambda_BA.4_5, lambda_BQ, lambda_XBB, lambda_EG_HV, lambda_JN, lambda_KP_LB))
}

run_metropolis_MCMC <- function(startvalue, iterations){
  chain = array(dim = c(iterations+1,20))
  chain[1,] = c(startvalue, posterior(proposalfunction(startvalue)), 0)

    for (i in 1:iterations){
      
    proposal = proposalfunction(chain[i,1:18])
    
    print("Proposal: ")
    print(proposal)
    
    # old <- Sys.time()
    
    posterior_proposal <- posterior(proposal)
    # new <- Sys.time() - old # calculate difference
    # print(new)
  
    
    posterior_current <- chain[i,19]
    # print(paste0("Posterior proposal: ", posterior_proposal))
    # print(paste0("Posterior current: ", posterior_current))
    
    acceptance_ratio <- sigmoid((posterior_proposal - posterior_current))
    if(posterior_proposal == Inf){
      acceptance_ratio <- 0
    }
    random_probab <- runif(1)
    
    # print(paste0("Probability: ", acceptance_ratio))
    # print(paste0("Random prob: ", random_probab))
    
    #Using sigmoid function for acceptance ratio
    if (random_probab < acceptance_ratio){
      chain[i+1,1:18] <- proposal
      chain[i+1,19] <- posterior_proposal
      chain[i+1,20] <- acceptance_ratio
      
    }else{
      chain[i+1,1:18] <- chain[i, 1:18]
      chain[i+1,19] <- posterior_current
      chain[i+1,20] <- acceptance_ratio
    }
    
    # print(paste0("Iteration ", i, " Chain: "))
    # print(chain[i+1,])
  }
  return(chain)
}

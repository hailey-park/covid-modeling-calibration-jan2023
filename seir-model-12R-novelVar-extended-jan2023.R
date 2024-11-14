###################################################################################################
#Title: SEIR Model Function
#Author: Hailey Park
#Date: May 23, 2024
###################################################################################################

#Read in data
beta_coeffs <- read.csv("data/clean-data/beta_coeffs.csv")[,-1]   #Beta coefficients (age-specific, immunity type-specific, and waning group-specific)    
severity_coeffs <- read.csv("data/clean-data/severity_coeffs.csv")[,-1]  #Severity coefficients (age-specific, immunity type-specific, and waning group-specific)    
primary_vax_rate <- read.csv("data/clean-data/primary_vax_rate_over_time_extended_oct2024.csv")[,-1]  #Primary series rate over time by age group
booster_vax_rate <- read.csv("data/clean-data/booster_vax_rate_over_time_extended_oct2024.csv")[,-1]  #Booster rate over time by age group
contact_matrix <- read.csv("data/clean-data/contact matrix.csv")
novel_variants_share <- read.csv("data/clean-data/variant_data_clean_extended_oct2024.csv")[,-1]
lambda_df <- data.frame(variant_group = c("BA.1.1", "BA.2", "BA.4/5", "BQ", "XBB", "EG/HV", "JN.1", "KP/LB"),
                        lambda = rep(0, 8))


#Adjusting the time range for some data (April 2022 start date instead of January 2022)
primary_vax_rate <- primary_vax_rate %>% filter(week >= as.Date("2022-12-25")) %>% mutate(week_num = week_num - 52)
booster_vax_rate <- booster_vax_rate %>% filter(week >= as.Date("2022-12-25")) %>% mutate(week_num = week_num - 52)
novel_variants_share <- novel_variants_share %>% filter(week >= as.Date("2022-12-24")) %>% mutate(weeks_since = weeks_since - 52)

# Function for running age-structured SEIR model of COVID transmission
# Features: waning immunity, different immunity types
# Age Groups: 0-17, 18-49, 50-64, 65-74, 75+ years
seir_model = function (current_timepoint, state_values, parameters)
{
  contact_matrix_adj <- c(parameters[19:23])
  
  parameters <- c(severity_0_17 = parameters[1],
                  severity_18_49 = parameters[2],
                  severity_50_64 = parameters[3],
                  severity_65_74 = parameters[4],
                  severity_75plus = parameters[5],
                  nonsevere_waning_rate_vaccine = parameters[6],
                  nonsevere_waning_rate_hybrid = parameters[7],
                  lambda = parameters[8], #this is ignored for now
                  severe_waning_level_vaccine = parameters[9],
                  severe_waning_level_hybrid = parameters[10],
                  lambda_BA.1.1 = parameters[11],
                  lambda_BA.2 = parameters[12],
                  lambda_BA.4_5 = parameters[13],
                  lambda_BQ = parameters[14],
                  lambda_XBB = parameters[15],
                  lambda_EG_HV = parameters[16],
                  lambda_JN = parameters[17],
                  lambda_KP_LB = parameters[18])
  n_ages <- 5

  # create state variables (local variables)
  #Immune naive
  S_0 <- as.matrix(state_values[1:n_ages])                 # susceptible
  E_0 <- as.matrix(state_values[(n_ages+1):(2*n_ages)])    # exposed
  I_0 <- as.matrix(state_values[(2*n_ages+1):(3*n_ages)])  # infectious
  H_0 <- as.matrix(state_values[(3*n_ages+1):(4*n_ages)])  # hospitalization/death
  
  #Vaccine-induced immunity only (with waning)
  S_V_1 <- as.matrix(state_values[(4*n_ages+1):(5*n_ages)])     #waning 0-3 months
  E_V_1 <- as.matrix(state_values[(5*n_ages+1):(6*n_ages)]) 
  I_V_1 <- as.matrix(state_values[(6*n_ages+1):(7*n_ages)])
  H_V_1 <- as.matrix(state_values[(7*n_ages+1):(8*n_ages)]) 
  
  S_V_2 <- as.matrix(state_values[(8*n_ages+1):(9*n_ages)])     #waning 3-6 months
  E_V_2 <- as.matrix(state_values[(9*n_ages+1):(10*n_ages)])
  I_V_2 <- as.matrix(state_values[(10*n_ages+1):(11*n_ages)])
  H_V_2 <- as.matrix(state_values[(11*n_ages+1):(12*n_ages)])
  
  S_V_3 <- as.matrix(state_values[(12*n_ages+1):(13*n_ages)])   #waning 6-12 months
  E_V_3 <- as.matrix(state_values[(13*n_ages+1):(14*n_ages)])
  I_V_3 <- as.matrix(state_values[(14*n_ages+1):(15*n_ages)])
  H_V_3 <- as.matrix(state_values[(15*n_ages+1):(16*n_ages)])
  
  S_V_4 <- as.matrix(state_values[(16*n_ages+1):(17*n_ages)])   #waning 12+ months
  E_V_4 <- as.matrix(state_values[(17*n_ages+1):(18*n_ages)])
  I_V_4 <- as.matrix(state_values[(18*n_ages+1):(19*n_ages)])
  H_V_4 <- as.matrix(state_values[(19*n_ages+1):(20*n_ages)])
  
  #Hybrid/Infection-induced immunity (with waning)
  S_1 <- as.matrix(state_values[(20*n_ages+1):(21*n_ages)])     #waning 0-3 months
  E_1 <- as.matrix(state_values[(21*n_ages+1):(22*n_ages)])
  I_1 <- as.matrix(state_values[(22*n_ages+1):(23*n_ages)])
  H_1 <- as.matrix(state_values[(23*n_ages+1):(24*n_ages)])

  S_2 <- as.matrix(state_values[(24*n_ages+1):(25*n_ages)])     #waning 3-6 months
  E_2 <- as.matrix(state_values[(25*n_ages+1):(26*n_ages)])
  I_2 <- as.matrix(state_values[(26*n_ages+1):(27*n_ages)])
  H_2 <- as.matrix(state_values[(27*n_ages+1):(28*n_ages)]) 
  
  S_3 <- as.matrix(state_values[(28*n_ages+1):(29*n_ages)])     #waning 6-12 months
  E_3 <- as.matrix(state_values[(29*n_ages+1):(30*n_ages)])
  I_3 <- as.matrix(state_values[(30*n_ages+1):(31*n_ages)])
  H_3 <- as.matrix(state_values[(31*n_ages+1):(32*n_ages)])
  
  S_4 <- as.matrix(state_values[(32*n_ages+1):(33*n_ages)])     #waning 12+ months
  E_4 <- as.matrix(state_values[(33*n_ages+1):(34*n_ages)])
  I_4 <- as.matrix(state_values[(34*n_ages+1):(35*n_ages)])
  H_4 <- as.matrix(state_values[(35*n_ages+1):(36*n_ages)])
  
  #Recovered compartment -- for perfect immunity after infection
  R_1 <- as.matrix(state_values[(36*n_ages+1):(37*n_ages)])       # recovered stage 1 
  R_2 <- as.matrix(state_values[(37*n_ages+1):(38*n_ages)])       # recovered stage 2 
  R_3 <- as.matrix(state_values[(38*n_ages+1):(39*n_ages)])       # recovered stage 3
  R_4 <- as.matrix(state_values[(39*n_ages+1):(40*n_ages)])       # recovered stage 4
  R_5 <- as.matrix(state_values[(40*n_ages+1):(41*n_ages)])       # recovered stage 5
  R_6 <- as.matrix(state_values[(41*n_ages+1):(42*n_ages)])       # recovered stage 6
  R_7 <- as.matrix(state_values[(42*n_ages+1):(43*n_ages)])       # recovered stage 7
  R_8 <- as.matrix(state_values[(43*n_ages+1):(44*n_ages)])       # recovered stage 8
  R_9 <- as.matrix(state_values[(44*n_ages+1):(45*n_ages)])       # recovered stage 9
  R_10 <- as.matrix(state_values[(45*n_ages+1):(46*n_ages)])       # recovered stage 10
  R_11 <- as.matrix(state_values[(46*n_ages+1):(47*n_ages)])       # recovered stage 11
  R_12 <- as.matrix(state_values[(47*n_ages+1):(48*n_ages)])       # recovered stage 12 
    
  with ( 
    as.list (parameters),     # variable names within parameters can be used 
    {
      #total infections include nonsevere and severe infections
      total_inf <- I_0 + I_V_1 + I_V_2 + I_V_3 + I_V_4 + I_1 + I_2 + I_3 + I_4 + H_0 + H_V_1 + H_V_2 + H_V_3 + H_V_4 + H_1 + H_2 + H_3 + H_4
      
      #vaccination rates are time-varying to match CA vaccine administration
      current_primary_vax_rate <- primary_vax_rate[primary_vax_rate$week_num == round((current_timepoint/7)), ]$primary_series_rate
      current_booster_vax_rate <- booster_vax_rate[booster_vax_rate$week_num == round((current_timepoint/7) ), ]$boosted_rate

      if(is.na(current_primary_vax_rate[1])){
        current_primary_vax_rate <- rep(0, 5)
      }

      if(is.na(current_booster_vax_rate[1])){
        current_booster_vax_rate <- rep(0, 5)
      }
      
      #calculate current novel variant-adjusted lambda
      # novel_variants_with_lambdas <- merge(novel_variants_share %>% filter(weeks_since == max(ceiling((current_timepoint)/7), 1)),
      #                                      data.frame(variant_group = c("BA.1.1", "BA.2", "BA.4/5", "BQ", "XBB", "EG/HV", "JN.1", "KP/LB"),
      #                                                 lambda = c(lambda_BA.1.1, lambda_BA.2, lambda_BA.4_5, lambda_BQ, lambda_XBB, lambda_EG_HV, lambda_JN, lambda_KP_LB)),
      #                                      by = "variant_group", all.x = TRUE)
      # 
      #adj_lambda <- weighted.mean(novel_variants_with_lambdas$lambda, novel_variants_with_lambdas$rescaled_share)
      
      DF1 <- novel_variants_share %>% filter(weeks_since == max(ceiling((current_timepoint)/7), 1))
      # DF2 <- data.frame(variant_group = c("BA.1.1", "BA.2", "BA.4/5", "BQ", "XBB", "EG/HV", "JN.1", "KP/LB"),
      #                   lambda = c(lambda_BA.1.1, lambda_BA.2, lambda_BA.4_5, lambda_BQ, lambda_XBB, lambda_EG_HV, lambda_JN, lambda_KP_LB))
      # 
      
      lambda_df$lambda <- c(lambda_BA.1.1, lambda_BA.2, lambda_BA.4_5, lambda_BQ, lambda_XBB, lambda_EG_HV, lambda_JN, lambda_KP_LB)
      
      DF1$lambda = lambda_df$lambda[match(DF1$variant_group, lambda_df$variant_group)]
      adj_lambda <- weighted.mean(DF1$lambda, DF1$rescaled_share)
      

      #transmission coefficients by waning compartment and immunity type (fitted to nonsevere protection curves)
      # age-stratified (5 age groups)
      # calculation: beta coefficient = (1 - average nonsevere protection by waning compartment + immunity type)
      beta_0 <- 1              
      beta_V_1 <- beta_coeffs$beta_V_1
      beta_V_2 <- beta_coeffs$beta_V_2
      beta_V_3 <- beta_coeffs$beta_V_3
      beta_V_4 <- beta_coeffs$beta_V_4
      beta_1 <- beta_coeffs$beta_1
      beta_2 <- beta_coeffs$beta_2
      beta_3 <- beta_coeffs$beta_3
      beta_4 <- beta_coeffs$beta_4
      
      #severity case rates by waning compartment and immunity type (fitted to severe protection curves)
      # age-stratified (5 age groups)
      # severity case rates includes adjustment for nonsevere protection (beta term) because H compartment branch from E compartment
      # calculation: severity_rate = (1- average severe protection by waning compartment + immunity type)/(beta) * baseline severity rate
      severity_0 <- c(severity_0_17, severity_18_49, severity_50_64, severity_65_74, severity_75plus)       
      severity_V_1 <- ((severity_coeffs$severity_V_1 *  beta_V_1) + severe_waning_level_vaccine)/beta_V_1 * severity_0 
      severity_V_2 <- ((severity_coeffs$severity_V_2 *  beta_V_2) + severe_waning_level_vaccine)/beta_V_2 * severity_0 
      severity_V_3 <- ((severity_coeffs$severity_V_3 *  beta_V_3) + severe_waning_level_vaccine)/beta_V_3 * severity_0 
      severity_V_4 <- ((severity_coeffs$severity_V_4 *  beta_V_4) + severe_waning_level_vaccine)/beta_V_4 * severity_0 
      severity_1 <- ((severity_coeffs$severity_1 *  beta_1) + severe_waning_level_hybrid)/beta_1 * severity_0
      severity_2 <- ((severity_coeffs$severity_2 *  beta_2) + severe_waning_level_hybrid)/beta_2 * severity_0
      severity_3 <- ((severity_coeffs$severity_3 *  beta_3) + severe_waning_level_hybrid)/beta_3 * severity_0
      severity_4 <- ((severity_coeffs$severity_4 *  beta_4) + severe_waning_level_hybrid)/beta_4 * severity_0
      
      #waning rates by waning compartment and immunity type
      waning_rate_V_1 <- 1/90 * nonsevere_waning_rate_vaccine 
      waning_rate_V_2 <- 1/90 * nonsevere_waning_rate_vaccine 
      waning_rate_V_3 <- 1/180 * nonsevere_waning_rate_vaccine 
      waning_rate_1 <- 1/90 * nonsevere_waning_rate_hybrid 
      waning_rate_2 <- 1/90 * nonsevere_waning_rate_hybrid
      waning_rate_3 <- 1/180 * nonsevere_waning_rate_hybrid 
      waning_rate_R <- 1/(90/12) 
      
      #population count by age group
      pop_count <- total_pop_by_age$total_pop
      
      #contact matrix
      contact_matrix_clean <- as.matrix(contact_matrix[,-1] %>% t())
      dynamic_term <- (contact_matrix_clean %*% as.matrix(total_inf/pop_count))
      
      # compute derivatives
      dS_0 <- (-beta_0 * adj_lambda *  S_0 * dynamic_term * contact_matrix_adj) - (current_primary_vax_rate * S_0)
      dE_0 <- (beta_0 * adj_lambda * S_0 * dynamic_term * contact_matrix_adj) - (delta * E_0)
      dI_0 <- (delta * E_0 * (1 - severity_0)) - (gamma * I_0)
      dH_0 <- (delta * E_0 * severity_0) - (gamma * H_0)
      
      dS_V_1 <- (current_primary_vax_rate * S_0) + (current_booster_vax_rate * (S_V_2 + S_V_3 + S_V_4) - (beta_V_1 * adj_lambda * S_V_1 *  dynamic_term * contact_matrix_adj) - (waning_rate_V_1 * S_V_1))
      dE_V_1 <- (beta_V_1 * adj_lambda * S_V_1 * dynamic_term * contact_matrix_adj) - (delta * E_V_1)
      dI_V_1 <-  (delta * E_V_1 * (1 - severity_V_1)) - (gamma * I_V_1)
      dH_V_1 <- (delta * E_V_1 * severity_V_1) - (gamma * H_V_1)
      
      dS_V_2 <- (waning_rate_V_1 * S_V_1)  - (current_booster_vax_rate * S_V_2) - (waning_rate_V_2 * S_V_2) - (beta_V_2 * adj_lambda * S_V_2 * dynamic_term * contact_matrix_adj)
      dE_V_2 <- (beta_V_2 * adj_lambda * S_V_2 * dynamic_term * contact_matrix_adj) - (delta * E_V_2)
      dI_V_2 <- (delta * E_V_2 * (1 - severity_V_2)) - (gamma * I_V_2) 
      dH_V_2 <- (delta * E_V_2 * severity_V_2) - (gamma * H_V_2)

      dS_V_3 <- (waning_rate_V_2 * S_V_2) - (current_booster_vax_rate * S_V_3) - (waning_rate_V_3 * S_V_3) - (beta_V_3 * adj_lambda * S_V_3 * dynamic_term * contact_matrix_adj)
      dE_V_3 <- (beta_V_3 * adj_lambda * S_V_3 * dynamic_term * contact_matrix_adj) - (delta * E_V_3)
      dI_V_3 <- (delta * E_V_3 * (1 - severity_V_3)) - (gamma * I_V_3) 
      dH_V_3 <- (delta * E_V_3 * severity_V_3) - (gamma * H_V_3)
        
      dS_V_4 <- (waning_rate_V_3 * S_V_3) - (current_booster_vax_rate * S_V_4) - (beta_V_4 * adj_lambda * S_V_4 * dynamic_term * contact_matrix_adj)
      dE_V_4 <- (beta_V_4 * adj_lambda * S_V_4 * dynamic_term * contact_matrix_adj) - (delta * E_V_4)
      dI_V_4 <- (delta * E_V_4 * (1 - severity_V_4)) - (gamma * I_V_4)
      dH_V_4 <- (delta * E_V_4 * severity_V_4) - (gamma * H_V_4) 
        
      dS_1 <- (current_booster_vax_rate * (S_2 + S_3 + S_4)) - (beta_1 * adj_lambda * S_1 * dynamic_term * contact_matrix_adj) - (waning_rate_1 * S_1)
      dE_1 <- (beta_1 * adj_lambda * S_1 * dynamic_term * contact_matrix_adj) - (delta * E_1)
      dI_1 <- (delta * E_1 * (1 - severity_1)) - (gamma * I_1) 
      dH_1 <- (delta * E_1 * severity_1) - (gamma * H_1) 
        
      dS_2 <- (waning_rate_R * R_12) + (waning_rate_1 * S_1) - (current_booster_vax_rate * S_2) - (waning_rate_2 * S_2) - (beta_2 * adj_lambda * S_2 * dynamic_term * contact_matrix_adj)
      dE_2 <- (beta_2 * adj_lambda * S_2 * dynamic_term * contact_matrix_adj) - (delta * E_2)
      dI_2 <- (delta * E_2 * (1 - severity_2)) - (gamma * I_2) 
      dH_2 <- (delta * E_2 * severity_2) - (gamma * H_2) 
      
      dS_3 <- (waning_rate_2 * S_2) - (current_booster_vax_rate * S_3) - (waning_rate_3 * S_3) - (beta_3 * adj_lambda * S_3 * dynamic_term * contact_matrix_adj)
      dE_3 <- (beta_3 * adj_lambda * S_3 * dynamic_term * contact_matrix_adj) - (delta * E_3)
      dI_3 <- (delta * E_3 * (1 - severity_3)) - (gamma * I_3) 
      dH_3 <- (delta * E_3 * severity_3) - (gamma * H_3) 
        
      dS_4 <- (waning_rate_3 * S_3) - (current_booster_vax_rate * S_4) - (beta_4 * adj_lambda * S_4 * dynamic_term * contact_matrix_adj)
      dE_4 <- (beta_4 * adj_lambda * S_4 * dynamic_term * contact_matrix_adj) - (delta * E_4)
      dI_4 <- (delta * E_4 * (1 - severity_4)) - (gamma * I_4)
      dH_4 <- (delta * E_4 * severity_4) - (gamma * H_4)   
      
      dR_1 <- (gamma * total_inf) - (waning_rate_R * R_1)
      dR_2 <- (waning_rate_R * R_1) - (waning_rate_R * R_2)
      dR_3 <- (waning_rate_R * R_2) - (waning_rate_R * R_3)
      dR_4 <- (waning_rate_R * R_3) - (waning_rate_R * R_4)
      dR_5 <- (waning_rate_R * R_4) - (waning_rate_R * R_5)
      dR_6 <- (waning_rate_R * R_5) - (waning_rate_R * R_6)
      dR_7 <- (waning_rate_R * R_6) - (waning_rate_R * R_7)
      dR_8 <- (waning_rate_R * R_7) - (waning_rate_R * R_8)
      dR_9 <- (waning_rate_R * R_8) - (waning_rate_R * R_9)
      dR_10 <- (waning_rate_R * R_9) - (waning_rate_R * R_10)
      dR_11 <- (waning_rate_R * R_10) - (waning_rate_R * R_11)
      dR_12 <- (waning_rate_R * R_11) - (waning_rate_R * R_12)
      
      
      new_severe_cases <- (delta * E_0 * severity_0) + (delta * E_V_1 * severity_V_1) + (delta * E_V_2 * severity_V_2) + (delta * E_V_3 * severity_V_3) + (delta * E_V_4 * severity_V_4) + (delta * E_1 * severity_1) + (delta * E_2 * severity_2) + (delta * E_3 * severity_3) + (delta * E_4 * severity_4)
      #new_severe_cases <- (delta * E_0 * (1 - severity_0)) + (delta * E_V_1 * (1 - severity_V_1)) + (delta * E_V_2 * (1 - severity_V_2)) + (delta * E_V_3 * (1 - severity_V_3)) + (delta * E_V_4 * (1 - severity_V_4)) + (delta * E_1 * (1 - severity_1)) + (delta * E_2 * (1 - severity_2)) + (delta * E_3 * (1 - severity_3)) + (delta * E_4 * (1 - severity_4))

      # combine results
      results = c(dS_0, dE_0, dI_0, dH_0, 
                  dS_V_1, dE_V_1, dI_V_1, dH_V_1,
                  dS_V_2, dE_V_2, dI_V_2, dH_V_2,
                  dS_V_3, dE_V_3, dI_V_3, dH_V_3,
                  dS_V_4, dE_V_4, dI_V_4, dH_V_4,
                  dS_1, dE_1, dI_1, dH_1,
                  dS_2, dE_2, dI_2, dH_2,
                  dS_3, dE_3, dI_3, dH_3,
                  dS_4, dE_4, dI_4, dH_4,
                  dR_1, dR_2, dR_3, dR_4, dR_5, dR_6, dR_7, dR_8, dR_9, dR_10, dR_11, dR_12)
      list (results, new_severe_cases = new_severe_cases)
    }
  )
}






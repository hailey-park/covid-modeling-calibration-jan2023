###################################################################################################
#Title: Function for initializing population based on fitted severity parameters
#Author: Hailey Park
#Date: May 23, 2024
###################################################################################################

#Read in data
total_pop_by_age <- read.csv("data/clean-data/total_pop_by_age.csv")[,-1]   #Population count by age
beta_coeffs_relabeled <- read.csv("data/clean-data/beta_coeffs_relabeled.csv")[,-1]   #Beta coefficients (age-specific, immunity type-specific, and waning group-specific)    
severity_coeffs_relabeled <- read.csv("data/clean-data/severity_coeffs_relabeled.csv")[,-1]  #Severity coefficients (age-specific, immunity type-specific, and waning group-specific)    
severe_cases_total_init <- read.csv("data/clean-data/severe_cases_total_init_newDate_jan2023.csv")[,-1]  #Total severe cases in population (N = ~39 mil) by age at Jan 1 2022
contact_matrix <- read.csv("data/clean-data/contact matrix.csv")

# Function for initializing the population into correct compartments
model_pop_init = function (severity_params)
{
  # estimating proportion of severe cases out of all cases, with respect to each compartment (age group-, immunity type-, waning group-specific)
  # using baseline case hospitalization fraction in the immune naive population (severity_params)
  severity_proportion <- merge(severity_coeffs_relabeled, data.frame(age_group = c("0-17 years", "18-49 years", "50-64 years", "65-74 years", "75+ years"),
                                                           severity_param = severity_params), by = "age_group") %>% mutate(severity_prop = severity_coeff * severity_param)
  
  # merging the time-since results and beta terms to estimate relative proportions of susceptibles, exposed, and infected 
  dist_tbl_by_age <- merge(merge(merge(combined_time_since_by_age, beta_coeffs_relabeled, by = c("age_group", "immunity_type", "waning_group"), all.x = TRUE),
                           severity_proportion, by = c("age_group", "immunity_type", "waning_group"), all.x = TRUE),
                           total_pop_by_age, by = c("age_group")) %>% rename(relative_prop_susc = proportion, beta = value) %>%
    mutate(pop_count = total_pop * relative_prop_susc) %>% dplyr::select(c(age_group, immunity_type, waning_group, beta, severity_prop, pop_count)) %>%
     mutate(susc_X_beta = pop_count * beta,
            relative_prop_exposed = susc_X_beta/sum(susc_X_beta, na.rm = TRUE),
            exposed_X_severity = relative_prop_exposed * severity_prop,
            relative_prop_severe_inf = exposed_X_severity/sum(exposed_X_severity, na.rm = TRUE))
  
  # estimating total severe cases in each compartment
  abs_total_init_by_age <- merge(dist_tbl_by_age, severe_cases_total_init, by = "age_group") %>% group_by(age_group) %>%
    mutate(abs_total_severe_inf = relative_prop_severe_inf/sum(relative_prop_severe_inf, na.rm = TRUE) * severe_inf_total, 
           abs_total_nonsevere_inf = abs_total_severe_inf * (1-severity_prop)/severity_prop,
           abs_total_exposed = (abs_total_severe_inf + abs_total_nonsevere_inf)*3/5,
           abs_total_susc = pop_count - abs_total_severe_inf - abs_total_nonsevere_inf - abs_total_exposed,
           abs_total_recovered = if_else(immunity_type == "recovered", pop_count, 0)) %>%
    mutate_at(vars(abs_total_severe_inf, abs_total_nonsevere_inf, abs_total_exposed, abs_total_susc, abs_total_recovered, pop_count), list(~round(., 3))) %>%
    dplyr::select(age_group, immunity_type, waning_group, abs_total_severe_inf, abs_total_nonsevere_inf, abs_total_exposed, abs_total_susc, abs_total_recovered, pop_count) %>% arrange(age_group)
  
  # check that all compartment assignments add up to 39,512,223 individuals
  sum(colSums(abs_total_init_by_age %>% ungroup() %>% dplyr::select(-c(age_group, immunity_type, waning_group, pop_count)), na.rm = TRUE))
  
  # check that there are no negative totals (susceptible compartments)
  if(any(abs_total_init_by_age$abs_total_susc < 0, na.rm = TRUE)) {
    print("Some 'Susceptible' compartments are negative.")
    return(NA)
    }
  
  # set initial population conditions (as fractions of entire pop)
  S_0 <- (abs_total_init_by_age %>% filter(immunity_type == "immune-naive"))$abs_total_susc             # susceptible -- immune-naive
  E_0 <- (abs_total_init_by_age %>% filter(immunity_type == "immune-naive"))$abs_total_exposed          # exposed
  I_0 <- (abs_total_init_by_age %>% filter(immunity_type == "immune-naive"))$abs_total_nonsevere_inf    # infected
  H_0 <- (abs_total_init_by_age %>% filter(immunity_type == "immune-naive"))$abs_total_severe_inf       # hospitalized
  S_V_1 <- (abs_total_init_by_age %>% filter(immunity_type == "vaccine-only", waning_group == "0-3 months"))$abs_total_susc            # -- vaccine-only (0-3 months)
  E_V_1 <- (abs_total_init_by_age %>% filter(immunity_type == "vaccine-only", waning_group == "0-3 months"))$abs_total_exposed
  I_V_1 <- (abs_total_init_by_age %>% filter(immunity_type == "vaccine-only", waning_group == "0-3 months"))$abs_total_nonsevere_inf
  H_V_1 <- (abs_total_init_by_age %>% filter(immunity_type == "vaccine-only", waning_group == "0-3 months"))$abs_total_severe_inf
  S_V_2 <- (abs_total_init_by_age %>% filter(immunity_type == "vaccine-only", waning_group == "3-6 months"))$abs_total_susc            # -- vaccine-only (3-6 months)
  E_V_2 <- (abs_total_init_by_age %>% filter(immunity_type == "vaccine-only", waning_group == "3-6 months"))$abs_total_exposed
  I_V_2 <- (abs_total_init_by_age %>% filter(immunity_type == "vaccine-only", waning_group == "3-6 months"))$abs_total_nonsevere_inf
  H_V_2 <- (abs_total_init_by_age %>% filter(immunity_type == "vaccine-only", waning_group == "3-6 months"))$abs_total_severe_inf
  S_V_3 <- (abs_total_init_by_age %>% filter(immunity_type == "vaccine-only", waning_group == "6-12 months"))$abs_total_susc           # -- vaccine-only (6-12 months)
  E_V_3 <- (abs_total_init_by_age %>% filter(immunity_type == "vaccine-only", waning_group == "6-12 months"))$abs_total_exposed
  I_V_3 <- (abs_total_init_by_age %>% filter(immunity_type == "vaccine-only", waning_group == "6-12 months"))$abs_total_nonsevere_inf
  H_V_3 <- (abs_total_init_by_age %>% filter(immunity_type == "vaccine-only", waning_group == "6-12 months"))$abs_total_severe_inf
  S_V_4 <- rep(0, 5)  # -- vaccine-only (12+ months)
  E_V_4 <- rep(0, 5)
  I_V_4 <- rep(0, 5)
  H_V_4 <- rep(0, 5)
  S_1 <- (abs_total_init_by_age %>% filter(immunity_type == "hybrid/infection-only", waning_group == "0-3 months"))$abs_total_susc            # -- hybrid/infection-only (0-3 months)
  E_1 <- (abs_total_init_by_age %>% filter(immunity_type == "hybrid/infection-only", waning_group == "0-3 months"))$abs_total_exposed
  I_1 <- (abs_total_init_by_age %>% filter(immunity_type == "hybrid/infection-only", waning_group == "0-3 months"))$abs_total_nonsevere_inf
  H_1 <- (abs_total_init_by_age %>% filter(immunity_type == "hybrid/infection-only", waning_group == "0-3 months"))$abs_total_severe_inf
  S_2 <- (abs_total_init_by_age %>% filter(immunity_type == "hybrid/infection-only", waning_group == "3-6 months"))$abs_total_susc            # -- hybrid/infection-only (3-6 months)
  E_2 <- (abs_total_init_by_age %>% filter(immunity_type == "hybrid/infection-only", waning_group == "3-6 months"))$abs_total_exposed
  I_2 <- (abs_total_init_by_age %>% filter(immunity_type == "hybrid/infection-only", waning_group == "3-6 months"))$abs_total_nonsevere_inf
  H_2 <- (abs_total_init_by_age %>% filter(immunity_type == "hybrid/infection-only", waning_group == "3-6 months"))$abs_total_severe_inf
  S_3 <- (abs_total_init_by_age %>% filter(immunity_type == "hybrid/infection-only", waning_group == "6-12 months"))$abs_total_susc           # -- hybrid/infection-only (6-12 months)
  E_3 <- (abs_total_init_by_age %>% filter(immunity_type == "hybrid/infection-only", waning_group == "6-12 months"))$abs_total_exposed
  I_3 <- (abs_total_init_by_age %>% filter(immunity_type == "hybrid/infection-only", waning_group == "6-12 months"))$abs_total_nonsevere_inf
  H_3 <- (abs_total_init_by_age %>% filter(immunity_type == "hybrid/infection-only", waning_group == "6-12 months"))$abs_total_severe_inf
  S_4 <- (abs_total_init_by_age %>% filter(immunity_type == "hybrid/infection-only", waning_group == "12+ months"))$abs_total_susc             # -- hybrid/infection-only (12+ months)
  E_4 <- (abs_total_init_by_age %>% filter(immunity_type == "hybrid/infection-only", waning_group == "12+ months"))$abs_total_exposed
  I_4 <- (abs_total_init_by_age %>% filter(immunity_type == "hybrid/infection-only", waning_group == "12+ months"))$abs_total_nonsevere_inf
  H_4 <- (abs_total_init_by_age %>% filter(immunity_type == "hybrid/infection-only", waning_group == "12+ months"))$abs_total_severe_inf
  R_1 <- (abs_total_init_by_age %>% filter(immunity_type == "recovered"))$abs_total_recovered/12       # -- recovered
  R_2 <- (abs_total_init_by_age %>% filter(immunity_type == "recovered"))$abs_total_recovered/12       # -- recovered
  R_3 <- (abs_total_init_by_age %>% filter(immunity_type == "recovered"))$abs_total_recovered/12       # -- recovered
  R_4 <- (abs_total_init_by_age %>% filter(immunity_type == "recovered"))$abs_total_recovered/12       # -- recovered
  R_5 <- (abs_total_init_by_age %>% filter(immunity_type == "recovered"))$abs_total_recovered/12       # -- recovered
  R_6 <- (abs_total_init_by_age %>% filter(immunity_type == "recovered"))$abs_total_recovered/12       # -- recovered
  R_7 <- (abs_total_init_by_age %>% filter(immunity_type == "recovered"))$abs_total_recovered/12       # -- recovered
  R_8 <- (abs_total_init_by_age %>% filter(immunity_type == "recovered"))$abs_total_recovered/12       # -- recovered
  R_9 <- (abs_total_init_by_age %>% filter(immunity_type == "recovered"))$abs_total_recovered/12       # -- recovered
  R_10 <- (abs_total_init_by_age %>% filter(immunity_type == "recovered"))$abs_total_recovered/12       # -- recovered
  R_11 <- (abs_total_init_by_age %>% filter(immunity_type == "recovered"))$abs_total_recovered/12       # -- recovered
  R_12 <- (abs_total_init_by_age %>% filter(immunity_type == "recovered"))$abs_total_recovered/12       # -- recovered
  
  initial_values <- c(S_0 = S_0, E_0 = E_0, I_0 = I_0, H_0 = H_0, 
                      S_V_1 = S_V_1, E_V_1 = E_V_1, I_V_1 = I_V_1, H_V_1 = H_V_1,
                      S_V_2 = S_V_2, E_V_2 = E_V_2, I_V_2 = I_V_2, H_V_2 = H_V_2,
                      S_V_3 = S_V_3, E_V_3 = E_V_3, I_V_3 = I_V_3, H_V_3 = H_V_3,
                      S_V_4 = S_V_4, E_V_4 = E_V_4, I_V_4 = I_V_4, H_V_4 = H_V_4,
                      S_1 = S_1, E_1 = E_1, I_1 = I_1, H_1 = H_1,
                      S_2 = S_2, E_2 = E_2, I_2 = I_2, H_2 = H_2,
                      S_3 = S_3, E_3 = E_3, I_3 = I_3, H_3 = H_3,
                      S_4 = S_4, E_4 = E_4, I_4 = I_4, H_4 = H_4,
                      R_1 = R_1, R_2 = R_2, R_3 = R_3, R_4 = R_4, R_5 = R_5, R_6 = R_6,
                      R_7 = R_7, R_8 = R_8, R_9 = R_9, R_10 = R_10, R_11 = R_11, R_12 = R_12) 
  
  return(initial_values)
}

# Function for calculating the contact matrix adjustment factors at model initialization
contact_matrix_init = function (state_values)
{
  n_ages <- 5
  
  # create state variables (local variables)
  #Immune naive
  I_0 <- as.matrix(state_values[(2*n_ages+1):(3*n_ages)])  # infectious
  H_0 <- as.matrix(state_values[(3*n_ages+1):(4*n_ages)])  # hospitalization/death
  
  #Vaccine-induced immunity only (with waning)
  #waning 0-3 months
  I_V_1 <- as.matrix(state_values[(6*n_ages+1):(7*n_ages)])
  H_V_1 <- as.matrix(state_values[(7*n_ages+1):(8*n_ages)]) 
  
  #waning 3-6 months
  I_V_2 <- as.matrix(state_values[(10*n_ages+1):(11*n_ages)])
  H_V_2 <- as.matrix(state_values[(11*n_ages+1):(12*n_ages)])
  
  #waning 6-12 months
  I_V_3 <- as.matrix(state_values[(14*n_ages+1):(15*n_ages)])
  H_V_3 <- as.matrix(state_values[(15*n_ages+1):(16*n_ages)])
  
  #waning 12+ months
  I_V_4 <- as.matrix(state_values[(18*n_ages+1):(19*n_ages)])
  H_V_4 <- as.matrix(state_values[(19*n_ages+1):(20*n_ages)])
  
  #Hybrid/Infection-induced immunity (with waning)
  #waning 0-3 months
  I_1 <- as.matrix(state_values[(22*n_ages+1):(23*n_ages)])
  H_1 <- as.matrix(state_values[(23*n_ages+1):(24*n_ages)])
  
  #waning 3-6 months
  I_2 <- as.matrix(state_values[(26*n_ages+1):(27*n_ages)])
  H_2 <- as.matrix(state_values[(27*n_ages+1):(28*n_ages)]) 
  
  #waning 6-12 months
  I_3 <- as.matrix(state_values[(30*n_ages+1):(31*n_ages)])
  H_3 <- as.matrix(state_values[(31*n_ages+1):(32*n_ages)])
  
  #waning 12+ months
  I_4 <- as.matrix(state_values[(34*n_ages+1):(35*n_ages)])
  H_4 <- as.matrix(state_values[(35*n_ages+1):(36*n_ages)])

  #total infections include nonsevere and severe infections
  total_inf <- I_0 + I_V_1 + I_V_2 + I_V_3 + I_V_4 + I_1 + I_2 + I_3 + I_4 + H_0 + H_V_1 + H_V_2 + H_V_3 + H_V_4 + H_1 + H_2 + H_3 + H_4
  
  #population count by age group
  pop_count <- total_pop_by_age$total_pop
  
  #contact matrix
  contact_matrix_clean <- as.matrix(contact_matrix[,-1] %>% t())
  
  #Create contact matrix adjustment factors
  contact_matrix_adj <- (total_inf/pop_count)/(contact_matrix_clean %*% as.matrix(total_inf/pop_count))
  
  # print((total_inf/pop_count))
  # print((contact_matrix_clean %*% as.matrix(total_inf/pop_count)))
  
  return(contact_matrix_adj)
}
### Simulate and save artificial data ####
simulate_participants <- function(par_row) {
  gen_model <- par_row$gen_model
  parameters <- par_row[,-c(1:6)]
  dir.create(paste0("gen_model_", gen_model), showWarnings = FALSE)
  setwd(paste0("gen_model_", gen_model))
  parameters <- parameters[, apply(parameters, 2, function(x) !all(is.na(x)))]
  parameters$theta1 <- 1

  set.seed(par_row$Nsimu+par_row$ntrials+20*str_length(gen_model))
  sim_df <- simulateRTConf(paramDf = parameters, model = gen_model,
                           n = par_row$ntrials,  # number of trials per cond and stimulus
                           delta=0.01,        # discretization step for simulation (in sec.)
                           maxrt=30, simult_conf = TRUE)         # maximum decision time for each trial
  sim_df <- sim_df %>% filter(response !=0) # if maxrt is exceeded; response is set to 0

  if (gen_model=="dynaViTE") {
    thetas_lower <- c(t(quantile(subset(sim_df, response==-1)$conf,
                                 probs=cumsum(c(parameters[paste("prob_rating", 1:4, sep="")])))))
    thetas_upper <- c(t(quantile(subset(sim_df, response==1)$conf,
                                 probs=cumsum(c(parameters[paste("prob_rating", 1:4, sep="")])))))
    #parameters <- parameters[!grepl("rating", names(parameters))]
    parameters$theta1 <- NULL
    parameters[paste("thetaLower", 1:4, sep="")] <- thetas_lower
    parameters[paste("thetaUpper", 1:4, sep="")] <- thetas_upper
    sim_df[sim_df$response==1, "rating"] <- as.numeric(cut(sim_df[sim_df$response==1, "conf"],
                                                           breaks = c(-Inf, thetas_upper, Inf),
                                                           include.lowest = TRUE))
    sim_df[sim_df$response==-1, "rating"] <- as.numeric(cut(sim_df[sim_df$response==-1, "conf"],
                                                            breaks = c(-Inf, thetas_lower, Inf),
                                                            include.lowest = TRUE))
  } else {
    thetas_lower <- c(t(quantile(subset(sim_df, response==1)$conf,
                                 probs=cumsum(c(parameters[paste("prob_rating", 1:4, sep="")])))))
    thetas_upper <- c(t(quantile(subset(sim_df, response==2)$conf,
                                 probs=cumsum(c(parameters[paste("prob_rating", 1:4, sep="")])))))
    #parameters <- parameters[!grepl("rating", names(parameters))]
    parameters$theta1 <- NULL
    parameters[paste("thetaLower", 1:4, sep="")] <- thetas_lower
    parameters[paste("thetaUpper", 1:4, sep="")] <- thetas_upper
    sim_df[sim_df$response==2, "rating"] <- as.numeric(cut(sim_df[sim_df$response==2, "conf"],
                                                           breaks = c(-Inf, thetas_upper, Inf),
                                                           include.lowest = TRUE))
    sim_df[sim_df$response==1, "rating"] <- as.numeric(cut(sim_df[sim_df$response==1, "conf"],
                                                            breaks = c(-Inf, thetas_lower, Inf),
                                                            include.lowest = TRUE))
  }
  parameters$gen_model <- gen_model
  parameters$Nsimu <- par_row$Nsimu
  parameters$ntrials <- par_row$ntrials

  sim_df$Nsimu <- par_row$Nsimu

  save(parameters, sim_df, file=paste0("gen_model_", gen_model, "_simu_", par_row$Nsimu, ".RData"))
  setwd("..")
}

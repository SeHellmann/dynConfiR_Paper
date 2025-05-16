#### Fit                                    ####
# Define a function to be parallelized over for fitting
fun_sim_fit <- function(simu_row) {
  # Get information for currect simulation
  gen_model <- simu_row$gen_model
  fit_model <- simu_row$fit_model
#  fit_model <- simu_row$fit_model
  n_simu <- as.numeric(simu_row$Nsimu)
  
  #setwd(paste0("gen_model_", gen_model))
  load(paste0("../simulated_data/gen_model_", gen_model, "/gen_model_", gen_model, "_simu_", n_simu, ".RData"))
  
  dir.create(paste0("gen_model_", gen_model), showWarnings = FALSE)
  setwd(paste0("gen_model_", gen_model))

  parameter_names <- c("a", "v1", "v2", "v3", "v4", "v5", "sv", "z", "sz", "t0", "st0", 
                "tau", "w", "sigvis", "svis", "lambda", "b", "wint", "wrt", "wx",
                "thetaLower1", "thetaLower2","thetaLower3", "thetaLower4",
                "thetaUpper1", "thetaUpper2","thetaUpper3", "thetaUpper4")
  print(paste("Starting fitting:", n_simu, " for generative model:", gen_model, " and fitted model:", fit_model))
  
  # Fit models and combine output
  sim_df$participant <- sim_df$Nsimu
  t00 <- Sys.time()
  res_fit <- fitRTConf(sim_df, model=fit_model, restr_tau = "simult_conf", 
                       grid_search=TRUE,precision=3, 
                       nRatings=5, useparallel=FALSE, opts = list(nAttempts=4, nRestarts=2),
                       logging = TRUE)
  fitting_time <- as.numeric(difftime(Sys.time(), t00, units="min"))
  
  res_fit[,setdiff(parameter_names, names(res_fit))] <- NA
  res_fit$gen_model <- gen_model 
  res_fit$fittingmins <- fitting_time
  res_fit$fit_model <-  fit_model
  res_fit$Nsimu <- n_simu
  save(res_fit, file=paste0("fitted_result_genmodel_", gen_model, "_fitmodel_", fit_model,"_", n_simu, ".RData"))
  print(paste("Finished fitting:", n_simu, "gen model:", gen_model, " fitted model:", fit_model))
  setwd("..")
  return(res_fit)
}







#### Compute true negative log  likelihood of the true generative model
compute_true_negloglik <- function(simu_row) {
  # Get information for currect simulation
  gen_model <- simu_row$gen_model
  n_simu <- as.numeric(simu_row$Nsimu)
  
  #setwd(paste0("gen_model_", gen_model))
  load(paste0("../simulated_data/gen_model_", gen_model, "/gen_model_", gen_model, "_simu_", n_simu, ".RData"))
  
  
  if (!grepl("RM", gen_model)) {
    simu_row$true_negLogLik <- -LogLikWEV(sim_df, paramDf = parameters, 
                                 model = gen_model, simult_conf = TRUE)
  } else {
    simu_row$true_negLogLik <- -LogLikRM(sim_df, paramDf = parameters, 
                                model = gen_model, time_scaled=TRUE)
  }
  
  return(simu_row)
}


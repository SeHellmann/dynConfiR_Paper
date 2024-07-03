#### Fit                                    ####
# Define a function to be parallelized over for fitting
fun_sim_fit <- function(simu_row) {
  # Get information for currect simulation
  gen_model <- simu_row$gen_model
#  fit_model <- simu_row$fit_model
  n_simu <- as.numeric(simu_row$Nsimu)
  
  #setwd(paste0("gen_model_", gen_model))
  load(paste0("../simulated_data/gen_model_", gen_model, "/gen_model_", gen_model, "_simu_", n_simu, ".RData"))
  
  dir.create(paste0("gen_model_", gen_model), showWarnings = FALSE)
  setwd(paste0("gen_model_", gen_model))

  outnames <- c("a", "v1", "v2", "v3", "v4", "v5", "sv", "z", "sz", "t0", "st0", "tau", 
                "w", "sigvis", "svis", "lambda", "b", "wint", "wrt", "wx",
                "thetaLower1", "thetaLower2","thetaLower3", "thetaLower4",
                "thetaUpper1", "thetaUpper2","thetaUpper3", "thetaUpper4", 
                "true_negLogLik", "gen_model", "ntrials")
  print(paste("Starting fitting:", n_simu, " for generative model:", gen_model))
  
  out <- data.frame(matrix(NA, nrow=2, ncol=length(outnames)))
  colnames(out) <- outnames
  
  # Fit models and combine output
  sim_df$participant <- sim_df$Nsimu
  t00 <- Sys.time()
  res_fit <- fitRTConf(sim_df, model=gen_model, restr_tau = "simult_conf", 
                       nRatings=5, useparallel=FALSE, opts = list(nAttempts=4, nRestarts=4),
                       logging = TRUE)
  fitting_time <- as.numeric(difftime(Sys.time(), t00, units="min"))
  if (gen_model=="dynaViTE") {
    true_negLogLik <- -LogLikWEV(sim_df, paramDf = parameters, 
                                 model = gen_model, simult_conf = TRUE)
  } else {
    true_negLogLik <- -LogLikRM(sim_df, paramDf = parameters, 
                                 model = gen_model, time_scaled=TRUE)
  }

  
  out[1, colnames(res_fit)] <- res_fit
  out$model <- gen_model 
  out$fittingmins <- fitting_time
  parameters$model <- "TrueGen"
  parameters$negLogLik <- true_negLogLik
  out[2, colnames(parameters)] <- parameters
  out$Nsimu <- n_simu
  save(out, file=paste0("fitted_result_fitted_model_", gen_model, "_", n_simu, ".RData"))
  print(paste("Finished fitting:", n_simu, "gen model:", gen_model))
  setwd("..")
  return(out)
}

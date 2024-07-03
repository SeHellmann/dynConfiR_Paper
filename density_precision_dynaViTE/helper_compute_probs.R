compute_probsdynaViTE <- function(data, paramDf, precision=1e-5) {
  #### Get information from paramDf ####
  nConds <- 3
  nRatings <- 3
  V <- c(t(paramDf[,paste("v",1:(nConds), sep = "")]))
  SV <- rep(paramDf$sv, nConds)
  S <- rep(1, nConds)
  
  ## Recover confidence thresholds
  thetas_upper <- c(-1e+32, t(paramDf[,paste("thetaUpper",1:(nRatings-1), sep = "")]), 1e+32)
  thetas_lower <- c(-1e+32, t(paramDf[,paste("thetaLower",1:(nRatings-1), sep="")]), 1e+32)
  
  data <-data %>% mutate(vth1 = case_when(.data$response==1 ~ thetas_upper[.data$rating],
                                          .data$response==-1 ~ thetas_lower[(.data$rating)]),
                         vth2 = case_when(.data$response==1 ~ thetas_upper[(.data$rating+1)],
                                          .data$response==-1 ~ thetas_lower[(.data$rating+1)]),
                         M_drift = V[.data$condition]*.data$stimulus,
                         SV = SV[.data$condition],
                         S = S[.data$condition])
  
  
  probs <- with(data, dWEV(rt, vth1,vth2,
                           response=response,
                           tau=paramDf$tau, a=paramDf$a,
                           v = M_drift,
                           t0 = paramDf$t0, z = paramDf$z, sz = paramDf$sz, st0=paramDf$st0,
                           sv = SV, w=paramDf$w, muvis=abs(M_drift), svis=paramDf$svis,
                           sigvis=paramDf$sigvis, lambda=paramDf$lambda, s = S,
                           simult_conf = FALSE, z_absolute = FALSE,
                           precision = precision, stop_on_error = FALSE,
                           stop_on_zero=FALSE))
  return(probs)
}


compute_probsRM <-  function(data, paramDf, model="IRM", precision=NULL) {
  #### Get information from paramDf ####
  nConds <- 3
  nRatings <- 3
  V <- c(t(paramDf[,paste("v",1:(nConds), sep = "")]))
  
  thetas_1 <- c(1e-32, t(paramDf[,paste("thetaUpper",1:(nRatings-1), sep = "")]), 1e+64)
  thetas_2 <- c(1e-32, t(paramDf[,paste("thetaLower",1:(nRatings-1), sep="")]), 1e+64)
  
  ## Compute the row-wise likelihood of observations
  data <-data %>% mutate(a = paramDf$a,
                         b = paramDf$b,
                         th1 = case_when(.data$response==1 ~ thetas_1[(.data$rating)],
                                         .data$response==2 ~ thetas_2[(.data$rating)]),
                         th2 = case_when(.data$response==1 ~ thetas_1[(.data$rating+1)],
                                         .data$response==2 ~ thetas_2[(.data$rating+1)]),
                         mu1 = V[.data$condition]*(-1)^(1+.data$stimulus),
                         mu2 = V[.data$condition]*(-1)^(.data$stimulus),
                         t0 = paramDf$t0,
                         st0 = paramDf$st0)
  data$wx = paramDf$wx
  data$wrt = paramDf$wrt
  data$wint = paramDf$wint
  data$s = 1
  if (model=="IRMt") {
    probs <- dIRM(data$rt, data$response,data$mu1, data$mu2, data$a, data$b,
                  data$th1, data$th2, data$wx,  data$wrt,  data$wint,
                  data$t0, data$st0, data$s, step_width = precision)
  } else {
    probs <- dPCRM(data$rt, data$response,data$mu1, data$mu2, data$a, data$b,
                   data$th1, data$th2, data$wx,  data$wrt,  data$wint,
                   data$t0, data$st0, data$s, step_width = precision)
  }
  return(probs)
}
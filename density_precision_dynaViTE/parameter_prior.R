parameter_sampling <- function(n, model){
  if (model=="dynaViTE") {
    
    ## Typical values derived from three experiments 
    a <-  runif(n, 1.4,  4)
    v1 <- runif(n, 0.01, 1)
    dv1 <- runif(n, 0.1, 1.2)
    dv2 <- runif(n, 0.1, 1.2)
    v2 <- v1+dv1
    v3 <- v2+dv2
    sv <- runif(n, 0.04, 1.3)
    z <-  runif(n, 0.4,  0.6)
    sz <- runif(n, 0,    0.8)
    t0 <- runif(n, 0,  1.5)
    st0 <-runif(n, 0.15,1)
    tau <-runif(n, 0,  1.5)
    lambda <- runif(n, 0, 1.3)
    sigvis <- runif(n, 0.2, 2.6)
    svis  <- runif(n, 0.01, 1.1)
    w <- runif(n, .2, .9)
    
    sz <- (min(z, (1-z))*2)*sz
    
    ## Randomly choose overall proportion of confidence
    ## judgments
    prob_rating2 <- runif(n, 0.1, 0.5)
    temp <- runif(n, 0.2, 0.7)
    prob_rating1 <- temp * (1-prob_rating2)
    prob_rating3 <- (1-temp) * (1-prob_rating2)
    lambda <- runif(n, 0, 1.3)
    sigvis <- runif(n, 0, 9)
    svis  <- runif(n, 0, 3.65)
    w <- runif(n, .1, 1)
    
    parameters <- data.frame(a, v1, v2, v3, sv, z, sz, t0, st0, 
                             tau, w, lambda, sigvis, svis,
                             prob_rating1, prob_rating2, prob_rating3)
    parameters$theta1 <- 1
    parameters$N <- 1:n
    parlist <- split(parameters, seq(nrow(parameters)))
    helper_fct <- function(parameters) {
      sim_df <- simulateRTConf(paramDf = parameters, model = "dynaViTE",
                               n = 1000,  # number of trials per cond and stimulus  
                               delta=0.01,        # discretization step for simulation (in sec.) 
                               maxrt=30, simult_conf = FALSE)         # maximum decision time for each trial
      
      thetas_lower <- c(t(quantile(subset(sim_df, response==-1)$conf,
                                   probs=cumsum(c(parameters[paste("prob_rating", 1:2, sep="")])))))
      thetas_upper <- c(t(quantile(subset(sim_df, response==1)$conf,
                                   probs=cumsum(c(parameters[paste("prob_rating", 1:2, sep="")])))))
      #parameters <- parameters[!grepl("rating", names(parameters))]
      parameters$theta1 <- NULL
      parameters[paste("thetaLower", 1:2, sep="")] <- thetas_lower
      parameters[paste("thetaUpper", 1:2, sep="")] <- thetas_upper
      return(parameters)
    }
    parlist <- lapply(parlist, helper_fct)
    
    parameters <- do.call(rbind, parlist)
    return(parameters)
  } else {
    ## Typical values derived from three experiments 
    absum <- runif(n, 0.4, 2.6)
    abr <- runif(n, 0.4, 1/0.4)
    abind <- rbinom(n, 1, 0.5)
    arel = abr/(1+abr)
    brel = abr/(1+abr)
    arel[abind] = 1-arel[abind]
    brel[!abind] = 1-arel[!abind]
    b <- absum*brel
    a <- absum*arel
    
    v1 <- runif(n, 0.01, 0.5)
    dv2 <- runif(n, 0.2,  1)
    dv3 <- runif(n, 0.2,  1)
    v2 <- v1+dv2
    v3 <- v2+dv3
    t0 <- runif(n, 0,  1.9)
    st0 <-runif(n, 0.15,1)
    
    wint <- runif(n, .2, .95)
    wrtrel <- runif(n, 0, 1)
    wrt <- (1-wint)*wrtrel
    wx <- (1-wint)*(1-wrtrel)
    
    ## Randomly choose overall proportion of confidence
    ## judgments
    prob_rating2 <- runif(n, 0.1, 0.5)
    temp <- runif(n, 0.2, 0.7)
    prob_rating1 <- temp * (1-prob_rating2)
    prob_rating3 <- (1-temp) * (1-prob_rating2)
    
    
    parameters <- data.frame(a, b, v1, v2, v3, t0, st0, 
                             wx, wrt, wint,
                             prob_rating1, prob_rating2, prob_rating3)
    parameters$theta1 <- 1
    parameters$N <- 1:n
    parlist <- split(parameters, seq(nrow(parameters)))
    helper_fct <- function(parameters) {
      sim_df <- simulateRTConf(paramDf = parameters, model = model,
                               n = 1000,  # number of trials per cond and stimulus  
                               delta=0.01,        # discretization step for simulation (in sec.) 
                               maxrt=30, simult_conf = FALSE)         # maximum decision time for each trial
      
      thetas_lower <- c(t(quantile(subset(sim_df, response==1)$conf,
                                   probs=cumsum(c(parameters[paste("prob_rating", 1:2, sep="")])))))
      thetas_upper <- c(t(quantile(subset(sim_df, response==2)$conf,
                                   probs=cumsum(c(parameters[paste("prob_rating", 1:2, sep="")])))))
      #parameters <- parameters[!grepl("rating", names(parameters))]
      parameters$theta1 <- NULL
      parameters[paste("thetaLower", 1:2, sep="")] <- thetas_lower
      parameters[paste("thetaUpper", 1:2, sep="")] <- thetas_upper
      return(parameters)
    }
    parlist <- lapply(parlist, helper_fct)
    
    parameters <- do.call(rbind, parlist)
    return(parameters)
  }
}

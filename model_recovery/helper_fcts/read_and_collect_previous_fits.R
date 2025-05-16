## Define and call function to not mess with the user's environment
temp_fct <- function() {
  load("prevfits/parfits_fixed_and_free_inclDDMConf.RData")
  fits_Law <- subset(parfits, model %in% c("dynaViTE", "IRMt", "PCRMt", "2DSD"))
  load("prevfits/collected_fitsNpredicts_KonfMask.RData")
  fits_RM_KonfMask <- subset(fits_RMmodels, model %in% c("IRMt", "PCRMt"))
  load("prevfits/collected_fitsNpredicts_bothRDKExperiments.RData")
  fits_RM_RDK <- subset(fits_RMmodels, model %in% c("IRMt", "PCRMt"))
  load("prevfits/modelfits_threeexps.RData")
  fits_dynaViTE <- subset(fits, model %in% c("dynaViTE", "2DSD"))
  
  load("prevfits/fitsandpredictions_Shekar2021.RData")
  fits_RM_Shekhar <- subset(parfits, grepl(pattern="RMt", model)) %>%
    select(where(~!all(is.na(.x)))) %>%
    mutate(participant = participant + 300)
  
  rm(list=setdiff(ls(), c("fits_RM_RDK", "fits_RM_KonfMask", "fits_dynaViTE", "fits_Law", "fits_RM_Shekhar")))
  # table(fits_Law$model)
  # table(fits_RM_KonfMask$model)
  # table(fits_RM_RDK$model)
  # table(fits_RM_Shekhar$model)
  # table(fits_Law$model)
  
  fits_RM_Law <- fits_Law %>% subset(grepl(model, pattern="RM"))  %>%
    select(where(~!all(is.na(.x)))) %>%
    mutate(participant = participant + 400)
  fits_RM_Shekhar <- fits_RM_Shekhar %>% 
    mutate(v5=v3,
           v4=(v3+v2)/2,
           v3=v2,
           v2=(v1+v2)/2)
  fits_RM <- rbind(
    select(fits_RM_Law, -starts_with("theta"), -c(N, k, BIC, AICc, AIC, negLogLik, fixed)),
    select(mutate(fits_RM_KonfMask,participant = participant + 200), -starts_with("theta"), -c(N, k, BIC, AICc, AIC, negLogLik)),
    select(fits_RM_RDK, -starts_with("theta"), -c(N, k, BIC, AICc, AIC, negLogLik)),
    select(fits_RM_Shekhar, -starts_with("theta"), -c(N, k, BIC, AICc, AIC, negLogLik,fixed))
  )
  table(fits_RM$participant)
  
  
  fits_dynavite_Law <- fits_Law %>% subset(model%in% c("dynaViTE", "2DSD"))  %>%
    select(where(~!all(is.na(.x)))) %>%
    mutate(participant = participant + 400)
  fits_dynaViTE_Shekhar <- fits_dynaViTE %>% filter(grepl(pattern="Shekhar",experiment)) %>% 
    mutate(v5=v3,
           v4=(v3+v2)/2,
           v3=v2,
           v2=(v1+v2)/2) %>%
    select(where(~!all(is.na(.x))))
  
  fits_dynaViTE <- rbind(
    mutate(select(fits_dynavite_Law, -starts_with("theta"), -c(N, k, BIC, AICc, AIC, negLogLik, fixed)), experiment="Law & Lee (2021)"),
    select(filter(fits_dynaViTE,grepl(pattern="Hellmann", experiment)), -starts_with("theta"), -c(N, k, BIC, AICc, AIC, negLogLik, fixed)),
    select(fits_dynaViTE_Shekhar, -starts_with("theta"), -c(N, k, BIC, AICc, AIC, negLogLik,fixed))
  )
  table(fits_dynaViTE$participant)
  fits_RM <- left_join(fits_RM, distinct(fits_dynaViTE[,c("participant", "experiment")]))
  fits_IRMt <- subset(fits_RM, model=="IRMt")
  fits_PCRMt <- subset(fits_RM, model=="PCRMt")
  head(fits_RM)
  head(fits_dynaViTE)
  fits_2DSD <- subset(fits_dynaViTE, model=="2DSD")
  fits_dynaViTE <- subset(fits_dynaViTE, model=="dynaViTE")
  
  fits_IRMt <- fits_IRMt[,c("model", "participant","experiment", setdiff(names(fits_IRMt),c("model", "participant","experiment") ))]
  fits_PCRMt <- fits_PCRMt[,c("model", "participant","experiment", setdiff(names(fits_PCRMt),c("model", "participant","experiment") ))]
  fits_dynaViTE <- fits_dynaViTE[,c("model", "participant","experiment", setdiff(names(fits_dynaViTE),c("model", "participant","experiment") ))]
  fits_2DSD <- fits_2DSD[,c("model", "participant","experiment", setdiff(names(fits_2DSD),c("model", "participant","experiment") ))]
  save(fits_IRMt, fits_PCRMt, fits_dynaViTE, fits_2DSD, 
       file="helper_fcts/collected_fits_models.RData")
  return(list(dynaViTE = fits_dynaViTE, `2DSD`=fits_2DSD,
              IRMt = fits_IRMt, PCRMt=fits_PCRMt))
}
fits <- temp_fct()
rm(temp_fct)

###########################################################################
#####         Script for the parameter recovery analysis            #######
###########################################################################


#### NOTE: The parameter recovery involves simulating and fitting a lot of 
####        data and models. This takes a lot of time!
####        

# Sebastian Hellmann, 15.05.2024

###   Preamble and imports    ####
rm(list = ls())
# # insert path by hand:
# script_path <- paste0("C:/Users/...", # insert right path, here
#                       "/parrec")
# or use RStudio to find script file path
script_path <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(script_path)

library(tidyverse)   # for data processing
library(dynConfiR)   # for fitting and predicting the models
library(parallel)    # to run fitting in parallel over many cores

source("helper_fcts/read_and_collect_previous_fits.R")
## Source functions for generating random rating proportions (from a dirichlet)
source("helper_fcts/fun_sample_rating_props.R")
## Source functions for simulating and savin artificial data
source("helper_fcts/helper_simulate_artificialdata_fivesteps.R")
## Source fitting function wrapper
source("helper_fcts/fitting_function_parrec_fivesteps.R")
#=====================

## Sample model parameters ###
set.seed(20052023)
par_samples  <- expand.grid(gen_model=c("dynaViTE", "IRMt", "PCRMt"),
                            ntrials = c(50, 100, 200, 500, 1000)) %>%
  group_by(gen_model, ntrials) %>%
  do(fits[[.data$gen_model]][sample(1:nrow(fits[[.data$gen_model]]), size=100, replace=TRUE),]) %>% ungroup()
par_samples  <- cbind(par_samples, gen_rating_proportions(nrow(par_samples)))
par_samples$sz <- pmin(par_samples$sz, pmin(par_samples$z, 1-par_samples$z)/2 - 1e-3)

## To visualize the distribution of sampled parameters
# for (genmodel in names(fits)) {
#   plt_parameters <- ggplot(pivot_longer(select(filter(par_samples,model==genmodel),where(~!all(is.na(.x))), -participant, -experiment,-model),
#                       cols=-"gen_model", names_to = "parameter",
#                       values_to="prior_sample")) +
#     geom_histogram(aes(x=prior_sample, fill=gen_model), position = "dodge", bins=20)+
#     facet_wrap(~parameter, scales="free")+
#     ggtitle(genmodel)
#   show(plt_parameters)
#   
# }

# Cast data frame of parameters to a list, but also keep a data frame
par_samples <- cbind(Nsimu = 1:nrow(par_samples), par_samples)
planned_fit <- par_samples
par_samples <- split(par_samples, seq(nrow(par_samples)))

# Create a directory to save all results
dir.create("saved_details", showWarnings = FALSE)
setwd("saved_details")
save(par_samples, planned_fit, file="all_par_samples_list.RData")

## Generate random observations using the sampled parameter sets ###
dir.create("simulated_data", showWarnings = FALSE)
setwd("simulated_data")
lapply(par_samples, simulate_participants)
setwd("..")


### Fit the models to the simulated data #### 
dir.create("fit_results", showWarnings = FALSE)
setwd("fit_results")
n.cores <- parallel::detectCores() - 1
cl <- makeCluster(n.cores, "SOCK", outfile = "")
clusterEvalQ(cl, library(dynConfiR))
clusterEvalQ(cl, library(tidyverse))
clusterExport(cl, c("fun_sim_fit"))
recovery_collected_results <- parLapplyLB(cl, par_samples,fun_sim_fit)
stopCluster(cl)

setwd(script_path)




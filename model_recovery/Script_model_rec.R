###########################################################################
#####           Script for the model recovery analysis              #######
###########################################################################


#### NOTE: The model recovery involves simulating and fitting a lot of 
####        data and models. This takes a lot of time!
####        

# Sebastian Hellmann, 10.04.2025

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
library(ggh4x)
library(ggpubr)
windowsFonts(Times=windowsFont("Times"))
dir.create("figures", showWarnings = FALSE)
# Colorblind-friendly color palette:
model_levels <-  c(          "PCRMt", "IRMt",              "dynaViTE", "2DSD")
model_colors <- c(          "#E69F00","#009E73",            "#D55E00", "#CC79A7")
custom_theme <-       theme_bw() +
  theme(plot.margin = margin(0, 0, 0, 0, "cm"),
        text = element_text(size=9, family="Times"),
        axis.text = element_text(size=9, family="Times", color="black"))

#load("modelfits_threeexps.RData")
source("helper_fcts/read_and_collect_previous_fits.R")

## Source functions for generating random rating proportions (from a dirichlet)
source("helper_fcts/fun_sample_rating_props.R")
## Source functions for simulating and saving artificial data
source("helper_fcts/helper_simulate_artificialdata_fivesteps.R")

## Source fitting function wrapper
source("helper_fcts/fitting_function_parrec_fivesteps.R")


#=====================
if (!file.exists("saved_details/all_par_samples_list.RData")) {
  ## Sample model parameters ###
#  set.seed(23102024)
  set.seed(12122024)
  
  par_samples  <- expand.grid(gen_model=c("dynaViTE", "2DSD", "IRMt", "PCRMt"),
                              ntrials = 50) %>% #c(100, 200, 500)
    group_by(gen_model, ntrials) %>%
    do(fits[[.data$gen_model]][sample(1:nrow(fits[[.data$gen_model]]), size=50, replace=TRUE),]) %>% ungroup()
  par_samples  <- cbind(par_samples, gen_rating_proportions(nrow(par_samples)))
  # Bound sz parameter away from boundary
  par_samples$sz <- ifelse(!is.na(par_samples$sz), 
                           pmin(par_samples$sz, pmin(par_samples$z, 1-par_samples$z)*2 - 1e-2), NA)
  
  # # To visualize the distribution of sampled parameters
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
  planned_fit <- par_samples[,c("gen_model", "ntrials", "Nsimu")] %>%
    right_join(expand.grid(Nsimu=1:nrow(par_samples), fit_model = c("dynaViTE", "2DSD", "IRMt", "PCRMt")))
  planned_fit <- split(planned_fit, seq(nrow(planned_fit)))
  par_samples <- split(par_samples, seq(nrow(par_samples)))
  
  
  
  # Create a directory to save all results
  dir.create("saved_details", showWarnings = FALSE)
  setwd("saved_details")
  
  save(par_samples, planned_fit, file="all_par_samples_list.RData")
} else {
  setwd("saved_details")
  load("all_par_samples_list.RData")
}


dir.create("simulated_data", showWarnings = FALSE)
if (length(list.files("simulated_data/"))==0) {
  ## Generate random observations using the sampled parameter sets ###
  setwd("simulated_data")
  lapply(par_samples, simulate_participants)
  setwd("..")
}


### Fit the models to the simulated data #### 
if (!file.exists("fit_results/list_w_all_recoveries.RData")) {
  dir.create("fit_results", showWarnings = FALSE)
  setwd("fit_results")
  n.cores <- 80# parallel::detectCores() - 5
  cl <- makeCluster(n.cores, "SOCK", outfile = "")
  clusterEvalQ(cl, library(dynConfiR))
  clusterEvalQ(cl, library(tidyverse))
  clusterExport(cl, c("fun_sim_fit"))
  recovery_collected_results <- parLapplyLB(cl,planned_fit,fun_sim_fit)
  stopCluster(cl)
  save(recovery_collected_results, file="list_w_all_recoveries.RData")
} else {
  load("fit_results/list_w_all_recoveries.RData")
}

setwd(script_path)

df_recovery_results <- data.frame()
rel_names <- c("negLogLik", "N", "k", "BIC", "AICc", "AIC", "gen_model", "fittingmins", "fit_model", "Nsimu")
for (i in 1:length(recovery_collected_results)) {
  df_recovery_results <- rbind(df_recovery_results,
                               recovery_collected_results[[i]][,rel_names])
}

df_recovery_results$model_class <- ifelse(grepl("RM", df_recovery_results$fit_model),
                                          "Race Models", "Diffusion Models")



p_fittime <- ggplot(df_recovery_results, aes(x=fit_model, fill=gen_model, y=fittingmins/60))+
  geom_boxplot(position="dodge")+
  scale_fill_manual(name="Generative Model",
                    breaks = model_levels,values=model_colors )+
  scale_x_discrete(name="Fitted Model")+
  scale_y_continuous(name="\nFitting Time [hours]")+
  # , breaks=time_breaks, #c(0.1, 1, 10, 100, 1000),
  #                      labels = parse(text=paste("10^", log(time_breaks,base=10), sep="")))+
  facet_grid2(cols=vars(model_class), scales="free", independent="y")+
  custom_theme+
  # theme(strip.background = element_blank(),
  #       strip.text = element_blank())+
  theme(axis.title.x = element_text(margin=margin(t=6, b=0, r=0, l=0, unit = "pt")))#log2_trans())+
#facet_nested(cols=vars(as.factor(precision>3)), scales = "free", independent = "y")
p_fittime

ggsave("figures/Fittingtimes.tiff",
       width = 12, height=5, units="cm",dpi=600)
ggsave("figures/Fittingtimes.eps",
       width = 15, height=5, units="cm",dpi=600, device = cairo_ps)
# ggsave("../../Draft/figures/model_recovery/Fittingtimes.eps",
#        width = 15, height=5, units="cm",dpi=600, device = cairo_ps)


# Compute overall PEP and bootstrapped PEP
set.seed(2201)
bootstrapped_PEP_all <- data.frame()
BMS_recovery <- data.frame()
if (!file.exists("group_BMS_results.RData")) {
  for (i in model_levels) {
    temp <- df_recovery_results %>% filter(gen_model==i) %>%
      mutate(model=fit_model, subject=Nsimu %% 50 + 1,
             subject=as.numeric(as.factor(subject))) 
    
    bootstrapped_PEP <- lapply(1:1000, function(x) {
      temp %>% filter(subject %in% sample(1:max(temp$subject), size=10)) %>%
        group_BMS_fits() %>% `$`("model_weights")  %>%
        as.data.frame() %>%
        rownames_to_column("fit_model") %>%
        mutate(gen_model=i) %>% select(fit_model, pep, gen_model) %>%
        mutate(bootstrap_i = x)
    })
    bootstrapped_PEP <- do.call(rbind, bootstrapped_PEP)
    bootstrapped_PEP_all <- rbind(bootstrapped_PEP_all, bootstrapped_PEP)
    
    temp <- temp %>% group_BMS_fits() %>% `$`("model_weights")  %>%
      as.data.frame() %>%
      rownames_to_column("fit_model") %>%
      mutate(gen_model=i)
    BMS_recovery <- rbind(BMS_recovery, temp)
  }
  save(BMS_recovery, bootstrapped_PEP_all, file="group_BMS_results.RData")
} else {
  load("group_BMS_results.RData")
}

print(BMS_recovery)

p_bootstrap_PEP <- ggplot(bootstrapped_PEP_all, aes(x=gen_model, y=pep, fill=fit_model))+
  geom_boxplot(position="dodge", color="black")+
  xlab("Generative Model")+ ylab("Protected Exceedance Probability")+
  scale_fill_manual(name="Fitted Model",
                    breaks = model_levels,values=model_colors )+
  scale_y_continuous(breaks=c(0, 0.2, 0.4, 0.6, 0.8, 1))+
  custom_theme
p_bootstrap_PEP
# ggsave("figures/BootstrapPEP_Recovery.tiff",
#        width = 12, height=5, units="cm",dpi=600)
# ggsave("figures/BootstrapPEP_Recovery.eps",
#        width = 12, height=5, units="cm",dpi=600, device = cairo_ps)





sbj_weights <- data.frame()
for (i in model_levels) {
  temp <- df_recovery_results %>% filter(gen_model==i) %>%
      mutate(model=fit_model, subject=Nsimu %% 50 + 1) %>%
      subject_modelweights() %>%
      as.data.frame() %>%
      mutate(gen_model=i) %>%
      arrange(across(all_of(c(i, as.character(model_levels[model_levels!=i]))))) %>%
      mutate(plt_order = n():1)
    sbj_weights <- rbind(sbj_weights, temp)
}

sbj_weights <- sbj_weights %>%
  pivot_longer(all_of(model_levels), names_to = "fitted", values_to = "weight")
p_BIC_part<- ggplot(sbj_weights, aes(x=plt_order, y=weight, fill=fitted))+
  geom_bar(stat = "identity", color="black", linewidth=0.2, width=1)+
  scale_fill_manual(name="Fitted Model",
                    breaks = model_levels,values=model_colors )+
  facet_nested(rows=c(vars("Generating Model"),vars(gen_model)))+
  scale_y_continuous(name="BIC weight", breaks=c(0,0.5 ,1))+
  scale_x_continuous(name="Participant (reordered)", breaks=c(1, 10, 20, 30, 40, 50), expand = c(0.001, 0.001))+
  theme_minimal()+custom_theme+
  theme(plot.margin = margin(0.01, 0, 3, 0))
#p_BIC_part



ggpubr::ggarrange(p_BIC_part, p_bootstrap_PEP, 
                  nrow=2, heights = c(0.6, 0.4), 
                  common.legend = TRUE, legend = "bottom") +
 ggpubr::bgcolor("white")+
 ggpubr::border("white")
ggsave("figures/Model_recovery_results.tiff",
       width = 16, height=24, units="cm",dpi=600)
ggsave("figures/Model_recovery_results.eps",
       width = 15, height=16, units="cm",dpi=600, device = cairo_ps)
# ggsave("../../Draft/figures/model_recovery/Model_recovery_results.eps",
#        width = 15, height=16, units="cm",dpi=600, device = cairo_ps)






### Do model selection based on AIC
# Compute overall PEP and bootstrapped PEP
set.seed(10042025)
bootstrapped_PEP_all_AIC <- data.frame()
BMS_recovery_AIC <- data.frame()
if (!file.exists("group_BMS_results_AIC.RData")) {
  for (i in model_levels) {
    temp <- df_recovery_results %>% filter(gen_model==i) %>%
      mutate(model=fit_model, subject=Nsimu %% 50 + 1,
             subject=as.numeric(as.factor(subject))) 
    
    bootstrapped_PEP <- lapply(1:1000, function(x) {
      temp %>% filter(subject %in% sample(1:max(temp$subject), size=10)) %>%
        group_BMS_fits(measure="AIC") %>% `$`("model_weights")  %>%
        as.data.frame() %>%
        rownames_to_column("fit_model") %>%
        mutate(gen_model=i) %>% select(fit_model, pep, gen_model) %>%
        mutate(bootstrap_i = x)
    })
    bootstrapped_PEP <- do.call(rbind, bootstrapped_PEP)
    bootstrapped_PEP_all_AIC <- rbind(bootstrapped_PEP_all_AIC, bootstrapped_PEP)
    
    temp <- temp %>% group_BMS_fits(measure = "AIC") %>% `$`("model_weights")  %>%
      as.data.frame() %>%
      rownames_to_column("fit_model") %>%
      mutate(gen_model=i)
    BMS_recovery_AIC <- rbind(BMS_recovery_AIC, temp)
  }
  save(BMS_recovery_AIC, bootstrapped_PEP_all_AIC, file="group_BMS_results_AIC.RData")
} else {
  load("group_BMS_results_AIC.RData")
}

print(BMS_recovery_AIC)

p_bootstrap_PEP_AIC <- ggplot(bootstrapped_PEP_all_AIC, aes(x=gen_model, y=pep, fill=fit_model))+
  geom_boxplot(position="dodge", color="black", size=0.4)+
  xlab("Generative Model")+ ylab("Protected Exceedance Probability")+
  scale_fill_manual(name="Fitted Model",
                    breaks = model_levels,values=model_colors )+
  scale_y_continuous(breaks=c(0, 0.2, 0.4, 0.6, 0.8, 1))+
  custom_theme
p_bootstrap_PEP_AIC
ggsave("figures/BootstrapPEP_Recovery_AIC.tiff",
       width = 12, height=5, units="cm",dpi=600)
ggsave("figures/BootstrapPEP_Recovery_AIC.eps",
       width = 15, height=9, units="cm",dpi=600, device = cairo_ps)
ggsave("../../Draft/figures/appendix/BootstrapPEP_Recovery_AIC.eps",
       width = 15, height=9, units="cm",dpi=600, device = cairo_ps)




# ggpubr::ggarrange(p_fittime, p_BIC_part, p_bootstrap_PEP, 
#                   nrow=3, heights = c(0.3, 0.45, 0.25), 
#                   common.legend = TRUE, legend = "bottom") +
#   ggpubr::bgcolor("white")+
#   ggpubr::border("white")
# ggsave("figures/All_Model_recovery_results.tiff",
#        width = 16, height=21, units="cm",dpi=600)
# ggsave("figures/All_Model_recovery_results.eps",
#        width = 16, height=21, units="cm",dpi=600, device = cairo_ps)

#source("Unused_Visualizations.R")

###########################################################################
#####             Script to check the precision of the              #######
#####               approximated density in dynaViTE                #######
###########################################################################

#### NOTE: In this script we load previous results by default, if they are
####          available. You can rerun the analysis, but for us, it took
####          3.5 hours on 21 cores to calculate the results. 
####      

# Sebastian Hellmann, 03.06.2024

###   Preamble and imports    ####

##### Script to check the precision of the approximated densities

###   Preamble and imports    ####
rm(list = ls())
# # insert path by hand:
# script_path <- paste0("C:/Users/PPA859/Documents", # insert right path, here
#                       "/density_precision_dynaViTE")
# # or use RStudio to find script file path
script_path <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(script_path)
{
  library(dynConfiR)
  library(tidyverse)
  library(parallel)
  library(ggpubr)
}
windowsFonts(Times=windowsFont("Times New Roman"))

#=====================


if (!file.exists("saved_dynaViTE2.RDATA")) {
  ## Prior distributions
  source("parameter_prior.R")
  ## Source function to compute probabilities over simulated data
  source("helper_compute_probs.R")
  if (!file.exists("simulated_data_precision.RData")) {
    if (!file.exists("saved_parameter_samples.RData")) {
      set.seed(21112023)
      par_samples <- expand.grid(model=c("dynaViTE", "IRMt", "PCRMt")) %>%
        group_by(model) %>%
        do(parameter_sampling(50, .data$model)) %>% ungroup()
      save(par_samples, file="saved_parameter_samples.RData")
    } else {
      load("saved_parameter_samples.RData")
    }
    names(par_samples)
    par_samples <- split(par_samples, seq(nrow(par_samples)))
    simulate_drop0s <- function(paramDf) subset(simulateRTConf(paramDf, n=100), response!=0)
    simulated_data <- lapply(par_samples, simulate_drop0s)
    save(simulated_data, par_samples, file="simulated_data_precision.RData")
  } else {
    load("simulated_data_precision.RData")
  }
  
  par_samples_df <- do.call(rbind, par_samples)
  res <- rbind(expand.grid(listind = which(par_samples_df$model=="dynaViTE"),
                           precision = c(0.5, 1, 1.5, 2, 2.5, 3, 3.5, 4, 5)),
               expand.grid(listind = which(grepl(par_samples_df$model, pattern="RM")),
                           precision = c(1e-2, 1e-3, 1e-4, 1e-5, 1e-6, 1e-7, 1e-8, 1e-12)))
  
  compute_prec <- function(X) {
    t01 <- Sys.time()
    model <- par_samples[[X[[1]]]]$model
    simu <- subset(simulated_data[[X[[1]]]], response !=0)
    print(paste0("Start:", paste(X, collapse=", ")))
    if (model=="dynaViTE") {
      return(list(probs=compute_probsdynaViTE(simu, par_samples[[X[[1]]]],
                                              precision=X[[2]]), time=difftime(Sys.time(), t01, units = "min")))
    } else {
      return(list(probs=compute_probsRM(simu, par_samples[[X[[1]]]], model,
                                        precision=X[[2]]), time=difftime(Sys.time(), t01, units = "min")))
    }
  }
  parallel::detectCores()
  
  res_dynaViTE <- res[res[,1] %in% which(par_samples_df$model=="dynaViTE"),]
  res_dynaViTE <- split(res_dynaViTE, seq(nrow(res_dynaViTE)))
  Ncores <- 23
  cl <- makeCluster(type = "SOCK", Ncores)
  clusterExport(cl, c("compute_prec", "par_samples", "simulated_data"), envir = environment())
  clusterEvalQ(cl, source("helper_compute_probs.R"))
  clusterEvalQ(cl, library(dynConfiR))
  clusterEvalQ(cl, library(tidyverse))
  
  
  on.exit(try(stopCluster(cl), silent = TRUE))
  
  #prob_res <- clusterApplyLB(cl, res, fun = compute_prec)
  t00 <- Sys.time()
  prob_res <- parLapplyLB(cl, res_dynaViTE, compute_prec)
  stopCluster(cl)
  time_used <- difftime(Sys.time(), t00, units = "hour")
  save(prob_res, res_dynaViTE, time_used, file="saved_dynaViTE2.RDATA")
} else {
  load("saved_dynaViTE2.RDATA")
}
time_used


## Get vector of computation times
times <- do.call(rbind, prob_res)[,2]
times <- as.numeric(times)


#prob_res1 <- matrix(NA, nrow=length(prob_res), ncol=600)
#times <- matrix(NA, nrow=length(res_dynaViTE), ncol=3)
for (i in 1:length(prob_res)) {
  prob_res[[i]] <- c(prob_res[[i]][[1]], rep(NA, 600-length(prob_res[[i]][[1]])))
}
prob_res <- as.matrix(do.call(rbind, prob_res), nrow=length(res_dynaViTE), ncol=600)

res_dynaViTE <- do.call(rbind, res_dynaViTE)
precisions <- unique(res_dynaViTE$precision)
listinds <- unique(res_dynaViTE$listind)
res2 <- expand.grid(listinds=listinds, precision=precisions, precision_ref = precisions)
for (i in 1:nrow(res2)) {
  indprob_res <- which(res_dynaViTE$listind==res2$listinds[i] & res_dynaViTE$precision==res2$precision[i])
  indprob_res_reference <- which(res_dynaViTE$listind==res2$listinds[i] & res_dynaViTE$precision==res2$precision_ref[i])
  
  # indprob_res_reference <- which(res$listind==listinds[i] & res$precision==precisions[length(precisions)])
  error <- prob_res[indprob_res,] - prob_res[indprob_res_reference,]
  res2[i, "meanabserror"] <- mean(abs(error), na.rm = TRUE)
  res2[i, "meanerror"] <- mean((error), na.rm = TRUE)
  res2[i, "meansquarederror"] <- mean((error)^2, na.rm = TRUE)
  res2[i, "min"] <- times[i]
}
for ( i in 1:nrow(res2)) res2[i, "issubseq"] <- 
  res2$precision_ref[i] == c(precisions, 100)[which(precisions==res2$precision[i])+1]

###### Rescale the precision to the current scale 
### the precision argument in the dynConfiR package was adapted 14/03/24 
##  based on these results
res2$precision2 <- res2$precision
res2$precision <- res2$precision2+3.5
res2$precision_ref2 <-  res2$precision_ref
res2$precision_ref <-  res2$precision_ref2+3.5



custom_theme <-       theme_bw() +
  theme(plot.margin = margin(0, 0, 0, 0, "cm"),
        text = element_text(size=9, family="Times"),
        panel.grid.minor = element_blank(),  # switch off minor gridlines
        axis.text = element_text(size=9, family="Times", color="black"))

p1 <- ggplot(filter(res2, issubseq),
             aes(x=as.factor(precision), y=meanabserror))+
  geom_boxplot()+
  scale_y_continuous(name="Mean absolute difference between\nprecision and (precision+0.5)",
                     trans = "log", breaks=10^(-c(4, 5, 6, 7, 8)), 
                     labels = parse(text=paste("10^-", c(4, 5, 6, 7, 8), sep="")))+
  scale_x_discrete(name="", limits=factor(unique(res2$precision)))+guides(x="none")+
  custom_theme
p1
p2 <- ggplot(filter(res2, precision_ref==max(res2$precision) & precision != max(res2$precision)),
             aes(x=as.factor(precision), y=meanabserror))+
  geom_boxplot()+
  scale_y_continuous(name="Mean absolute difference between\nprecision and (precision=8.5)",
                     trans = "log", breaks=10^(-c(4, 5, 6, 7, 8)), 
                     labels = parse(text=paste("10^-", c(4, 5, 6, 7, 8), sep="")))+
  scale_x_discrete(name="", limits=factor(unique(res2$precision)))+guides(x="none")+
  custom_theme
p2
time_breaks <- c(0.1, 1, 10, 100, 1000)
p3 <- ggplot(subset(res2, !is.na(min)), aes(x=as.factor(precision), y=min*60))+
  geom_boxplot()+
  scale_x_discrete(name="Precision argument")+
  scale_y_continuous(name="\nComputation time [sec]", 
                     trans = "log", breaks=time_breaks, #c(0.1, 1, 10, 100, 1000),
                     labels = parse(text=paste("10^", log(time_breaks,base=10), sep="")))+
  custom_theme+
  theme(axis.title.x = element_text(margin=margin(t=10, b=0, r=0, l=0, unit = "pt")))
p3
ggarrange(p1, p2, p3, nrow=3, heights = c(0.32, 0.32, 0.36)) +
  annotate("rect", xmin=0.493, xmax=0.59, ymin=0.025, ymax=0.955, 
           fill=NA, color="darkred", linewidth=1)+
  annotate("text", x=mean(c(0.493, 0.59)), y=0.988, vjust=1, hjust=0.5,
                          color="darkred", label="Default",
           family="Times")
# ggsave("figures/precision_results_dynaViTE.jpg", height=17, width=17, units="cm",
#        dpi = 600)
# ## Only for Manuscript
# ggsave("../../Draft/figures/precision/precision_results_dynaViTE.eps", height=17, width=17, units="cm",
#        dpi = 600, dev=cairo_ps)
save(prob_res,  res_dynaViTE, res2, times, file="saved_dynaViTE_final.RDATA")


# ### Additional, more detailed plot
# ggplot(subset(res2, precision<precision_ref), aes(x=factor(precision), y=log(meanabserror)))+
#   geom_boxplot()+
#   facet_grid(rows = vars(precision_ref))#, scales = "free_y"

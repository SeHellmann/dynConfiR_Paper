###########################################################################
#####             Script to check the precision of the              #######
#####               approximated density in Race Models             #######
###########################################################################

# Sebastian Hellmann, 04.04.2025

###   Preamble and imports    ####

##### Script to check the precision of the approximated densities

###   Preamble and imports    ####
rm(list = ls())
# # insert path by hand:
# script_path <- paste0("C:/Users/PPA859/Documents", # insert right path, here
#                       "/density_precision_RM")
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
dir.create("figures", showWarnings = FALSE)

#=====================


if (!file.exists("saved_RMs2.RDATA")) {
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
  res <- rbind(expand.grid(listind = which(grepl(par_samples_df$model, pattern="RM")),
                           precision = c(seq(2, 7, by=0.5),9)))
  
  compute_prec <- function(X) {
    t01 <- Sys.time()
    model <- par_samples[[X[[1]]]]$model
    simu <- subset(simulated_data[[X[[1]]]], response !=0)
    print(paste0("Start:", paste(X, collapse=", ")))
    if (model=="dynaViTE") {
      probs=compute_probsdynaViTE(simu, par_samples[[X[[1]]]],
                                  precision=X[[2]])
    } else {
      probs=compute_probsRM(simu, par_samples[[X[[1]]]], model,
                            precision=X[[2]])
    }
    time=difftime(Sys.time(), t01, units = "min")
    precision = X[[2]]
    save(time, probs, simu, model, precision, file=paste0("autosave/listind_", X[[1]], "_precision_", X[[2]], ".RData"))
    return(list(probs=probs, time=time))
  }
  parallel::detectCores()
  res_RMs <- res[res[,1] %in% which(par_samples_df$model!="dynaViTE"),]
  res_RMs$model <- par_samples_df$model[res_RMs[,1]]
  res_RMs <- split(res_RMs, seq(nrow(res_RMs)))
  Ncores <- 15
  cl <- makeCluster(type = "SOCK", Ncores)
  clusterExport(cl, c("compute_prec", "par_samples", "simulated_data"), envir = environment())
  clusterEvalQ(cl, source("helper_compute_probs.R"))
  clusterEvalQ(cl, library(dynConfiR))
  clusterEvalQ(cl, library(tidyverse))
  
  
  on.exit(try(stopCluster(cl), silent = TRUE))
  dir.create("autosave", showWarnings = FALSE)
  #prob_res <- clusterApplyLB(cl, res, fun = compute_prec)
  t00 <- Sys.time()
  prob_res <- parLapplyLB(cl, res_RMs, compute_prec)
  stopCluster(cl)
  time_used <- difftime(Sys.time(), t00, units = "hour")
  
  save(prob_res, res_RMs, time_used, file="saved_RMs2.RDATA")
} else {
  load("saved_RMs2.RDATA")
}

time_used

## Get vector of computation times
times <- do.call(rbind, prob_res)[,2]
times <- as.numeric(times)
sum(times) / 60 / 60

#prob_res1 <- matrix(NA, nrow=length(prob_res), ncol=600)
#times <- matrix(NA, nrow=length(res_dynaViTE), ncol=3)
for (i in 1:length(prob_res)) {
  prob_res[[i]] <- c(prob_res[[i]][[1]], rep(NA, 600-length(prob_res[[i]][[1]])))
}
prob_res <- as.matrix(do.call(rbind, prob_res), nrow=length(res_dynaViTE), ncol=600)


res_RMs <- do.call(rbind, res_RMs)
precisions <- unique(res_RMs$precision)
listinds <- unique(res_RMs$listind)
res2 <- expand.grid(listinds=listinds, precision=precisions, precision_ref = precisions)
res2$model <- ifelse(res2$listinds < 101, "IRMt", "PCRMt")
for (i in 1:nrow(res2)) {
  indprob_res <- which(res_RMs$listind==res2$listinds[i] & res_RMs$precision==res2$precision[i])
  indprob_res_reference <- which(res_RMs$listind==res2$listinds[i] & res_RMs$precision==res2$precision_ref[i])

  # indprob_res_reference <- which(res$listind==listinds[i] & res$precision==precisions[length(precisions)])
  error <- prob_res[indprob_res,] - prob_res[indprob_res_reference,]
  res2[i, "meanabserror"] <- mean(abs(error), na.rm = TRUE)
  res2[i, "meanerror"] <- mean((error), na.rm = TRUE)
  res2[i, "meansquarederror"] <- mean((error)^2, na.rm = TRUE)
  res2[i, "min"] <- times[i]
}
for ( i in 1:nrow(res2)) res2[i, "issubseq"] <- 
  res2$precision_ref[i] == c(precisions, 100)[which(precisions==res2$precision[i])+1]

custom_theme <-       theme_bw() +
  theme(plot.margin = margin(0, 0, 0, 0, "cm"),
        text = element_text(size=9, family="Times"),
        axis.text = element_text(size=9, family="Times", color="black"))
#res2 <- subset(res2, precision %in% c(2,3,4,5,6,7,9))
p1 <- ggplot(filter(res2, issubseq & meanabserror > 0),
             aes(x=as.factor(precision), y=meanabserror))+
  geom_boxplot()+
  scale_y_continuous(name="Mean absolute difference between\nprecision and (precision+0.5)",
                     trans = "log", breaks=10^(-c(7, 5, 3, 1)), #c(0.1, 1, 10, 100, 1000),
                     labels = parse(text=paste("10^-", c(7, 5, 3, 1), sep="")),
                     expand = expansion(mult=0.05, add=c(0, 1)))+
  scale_x_discrete(name="", limits=factor(unique(res2$precision)))+guides(x="none")+
  facet_grid(cols=vars(model))+
  custom_theme
p1
p2 <- ggplot(filter(res2, precision_ref==max(res2$precision) & precision != max(res2$precision)),
             aes(x=as.factor(precision), y=meanabserror))+
  geom_boxplot()+
  scale_y_continuous(name="Mean absolute difference between\nprecision and (precision=9)",
                     trans = "log", breaks=10^(-c(7, 5, 3, 1)), #c(0.1, 1, 10, 100, 1000),
                     labels = parse(text=paste("10^-", c(7, 5, 3, 1), sep="")))+
  scale_x_discrete(name="", limits=factor(unique(res2$precision)))+guides(x="none")+
  # stat_summary(fun = mean, geom = "errorbar", aes(ymax = ..y.., ymin = ..y..),
  #              width = .75, linewidth=1, linetype="dashed")+
  facet_grid(cols=vars(model))+
  custom_theme+
  theme(strip.background = element_blank(),
        strip.text = element_blank())

p2
time_breaks <- 10^c(-1, 0)#c(0.01, 0.1, 1, 10, 100, 1000)
p3 <- ggplot(subset(res2, !is.na(min)), aes(x=as.factor(precision), y=min*60))+
  geom_boxplot()+
  scale_x_discrete(name="Precision argument")+
  scale_y_continuous(name="\nComputation time [sec]",
                     trans = "log", breaks=time_breaks, #c(0.1, 1, 10, 100, 1000),
                     labels = parse(text=paste("10^", log(time_breaks,base=10), sep="")))+
  facet_grid(cols=vars(model))+
  custom_theme+
  theme(strip.background = element_blank(),
        strip.text = element_blank())+
  theme(axis.title.x = element_text(margin=margin(t=6, b=0, r=0, l=0, unit = "pt")))
p3
text_height = 0.95
rect_height = 0.944
ggpubr::ggarrange(p1, p2, p3, nrow=3, heights = c(0.32, 0.32, 0.36)) +
  annotate("rect", 
           xmin= (5.5+(2.7-2)*11/3)/50, 
           xmax= (5.5+(3.22-2)*11/3)/50, 
           ymin=0.01, ymax=rect_height, 
           fill=NA, color="darkred", linewidth=1)+
  annotate("rect", 
           xmin= (5.5+(5.75-2)*11/3)/50, 
           xmax= (5.5+(6.25-2)*11/3)/50, 
           ymin=0.01, ymax=rect_height, 
           fill=NA, color="darkred", linewidth=1)+
  annotate("rect", 
           xmin= (28.7+(2.7-2)*11/3)/50, 
           xmax= (28.7+(3.2-2)*11/3)/50, 
           ymin=0.01, ymax=rect_height, 
           fill=NA, color="darkred", linewidth=1)+
  annotate("rect", 
           xmin= (28.7+(5.75-2)*11/3)/50, 
           xmax= (28.7+(6.24-2)*11/3)/50, 
           ymin=0.01, ymax=rect_height, 
           fill=NA, color="darkred", linewidth=1)+
  annotate("text", x=(28.7+(6-2)*11/3)/50, y=text_height, vjust=0, hjust=0.5,
           color="darkred", label="Default for density",
           family="Times")+
  annotate("text", x=(28.7+(3-2)*11/3)/50, y=text_height, vjust=0, hjust=0.5,
           color="darkred", label="Default for fitting",
           family="Times")+
  annotate("text", x=(5.5+(3-2)*11/3)/50, y=text_height, vjust=0, hjust=0.5,
           color="darkred", label="Default for fitting",
           family="Times")+
  annotate("text", x=(5.5+(6-2)*11/3)/50, y=text_height, vjust=0, hjust=0.5,
           color="darkred", label="Default for density",
           family="Times")

# ggsave("figures/precision_results_RM.jpg", height=17, width=17, units="cm",
#        dpi = 600)
# ggsave("figures/precision_results_RM.eps", height=17, width=17, units="cm",
#        dpi = 600, dev=cairo_ps)
# ggsave("../../Draft/figures/precision/precision_results_RM.eps", height=17,
#        width=17, units="cm",
#        dpi = 600, dev=cairo_ps)
save(prob_res, res_RMs, res2, times, file="saved_RM_final.RDATA")





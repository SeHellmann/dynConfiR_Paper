###########################################################################
#####  Visualizing the results of the parameter recovery analysis   #######
###########################################################################

# Sebastian Hellmann, 03.04.2025

###   Preamble and imports    ####
rm(list = ls())
# # insert path by hand:
# script_path <- paste0("C:/Users/...", # insert right path, here
#                       "/parrec")
# or use RStudio to find script file path
script_path <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(script_path)
dir.create("figures", showWarnings = FALSE)
windowsFonts(Times=windowsFont("Times New Roman"))

library(tidyverse)
library(ggh4x)       # for nested_wrap (creates better strip labels)

# Define parameter names and respective expressions for the labels
par_levels <- c(paste0("v", 1:5), 'sv', 'a', 'z', 'sz','b',
                't0', 'st0', 'tau', 'w', 'sigvis', 'svis',
                "lambda",
                'wx', 'wint', 'wrt',
                paste0('thetaLower',1:4), paste0("thetaUpper", 1:4))
par_labels <- c(paste0("nu[", 1:5, "]"), 's[nu]', 'a', 'z', 's[z]',
                'b',
                't[0]', 's[t0]', 'tau', 'w', 'sigma[Vis]', 's[Vis]', 
                'lambda', paste0("w[", c("x", "int", "rt"), "]"),
                paste0('theta[', -1,'~', 1:4,']'),#'',
                paste0('theta[', 1,'~',  1:4,']'))
parameter_classes <- data.frame(parameter=par_labels, 
                                class=factor(c(rep("Decision",12), rep("Confidence", 16)), levels=c("Decision", "Confidence"), ordered=TRUE))
models <- c("dynaViTE", "IRMt", "PCRMt")


load("saved_results.RData")
pars <- do.call(bind_rows, recovery_collected_results)

pars<- pars %>% 
  group_by(Nsimu) %>%
  mutate(ntrials = max(ntrials, na.rm=TRUE),
         model = model[1],
         class = ifelse(is.na(gen_model), "est", "true")) %>%
  ungroup() %>%
  filter(ntrials %in% c(50, 100, 200, 500))
#pars <- pars %>% mutate(wrt = ifelse(model=="IRMt" & class=="true", 1-wint-wx, wrt))

head(pars[,c("model", "fittingmins", "Nsimu", "class")])
table(pars$model, pars$ntrials, pars$class)

collected_par_recovery <- pars %>%
  select(-c(true_negLogLik, gen_model, fixed, N, k, BIC, AICc, AIC), -starts_with("prob_"))%>%
  pivot_longer(cols = a:thetaUpper4,
               names_to = "parameter", values_to = "value") %>%
  filter(!is.na(value)) %>%
  rename(n_sim=Nsimu) %>%
  arrange(n_sim, model, parameter, class) %>%
  pivot_wider(id_cols = c("n_sim", "parameter", "model", "ntrials"), names_from = class, values_from = value) %>%
  mutate(parameter = factor(parameter, levels =par_levels,labels = par_labels ))


## Compute the concordance correlation coefficient 
## between true and recovered parameters
CCC <- function(x, y) { ## Function to compute the 
  # concordance correlation coefficient 
  NAinds <- is.na(x) | is.na(y)
  x <- x[!NAinds]; y <- y[!NAinds]
  k <- length(y)
  sy2 <- var(y) * (k - 1)/k
  sx2 <- var(x) * (k - 1)/k
  est <- 2 * cov(x, y) *(k-1)/k /(sx2 + sy2 + (mean(y) - mean(x))^2)
  return(est)
}

est_cors <- collected_par_recovery %>%
  group_by(model, parameter, ntrials) %>%
  summarise(cor = CCC(est, true),
            x_min=min(true),
            y_max=max(est)) %>%
  left_join(parameter_classes) 
est_cors <- filter(est_cors, !is.na(cor))

filter(est_cors, parameter=="s[Vis]")

ggplot(est_cors, aes(x=parameter, y=cor, col=model))+ #shape=model, 
  geom_point(size=3)+
  facet_nested(cols=c(vars("Trials per condition and and stimulus"),vars(ntrials)),
               rows=c(vars("Model"), vars(model),vars(class)),scales = "free_x", independent = "x")+
  ylab("Concordance correlation coefficient")+
  scale_x_discrete(name="Parameter", labels = ggplot2:::parse_safe)+
  theme_bw()+ scale_color_discrete(guide="none")+
  ylim(c(min(est_cors$cor, na.rm=TRUE)-0.05, 1))+
  theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5))
ggsave(file = "figures/all_coeff.eps", width = 23.62, height=26, 
       units = "cm", dpi=600, device = cairo_ps)
# ## Only for manuscript
# ggsave(file = "../../Draft/figures/recovery/all_coeff.eps",
#        width = 23.62, height=26,
#        units = "cm", dpi=600, device = cairo_ps)


for (cur_model in unique(collected_par_recovery$model)) {
  for (cur_ntrials in unique(collected_par_recovery$ntrials)) {
    p_rec_plt <- ggplot(subset(collected_par_recovery, model==cur_model & ntrials==cur_ntrials), aes(x=true, y=est))+
      geom_smooth(data=subset(collected_par_recovery, model==cur_model&ntrials==cur_ntrials), method="lm", alpha=0.5, linewidth=0.9,
                  aes(x=true, y=est), formula = 'y~x', inherit.aes = FALSE)+
      geom_abline(col="black", linewidth=1, alpha=0.2)+
      geom_point(alpha=0.5)+
      scale_color_manual(breaks = c(TRUE, FALSE), name="Fitting stuck \nnear start value",
                         values=c("red", "black"))+
      geom_text(data = mutate(subset(est_cors, model==cur_model & ntrials == cur_ntrials)),
                aes(x=x_min, y=y_max,
                    label=paste0(format(round(cor, 2), nsmall=2))),#"rho == ",
                color="black",
                parse=TRUE, hjust=0, vjust=1)+
      facet_wrap(.~parameter, scales = "free", labeller = label_parsed,
                 drop = TRUE,  ncol=4)+
      theme_bw()+ylab("Fitted parameter value")+xlab("True parameter value")+
      ggtitle(paste0("Parameter recovery for ", cur_model," with ", cur_ntrials,"x 10 trials"))+
      geom_line(data=data.frame(parameter=rep(c("theta[-1~1]","theta[1~1]"),2),
                                est=c(1,1,1,1), true=c(0.6, 0.6, 1.4, 1.4)),
                col="blue")+
      theme(plot.margin = margin(0, 0.5, 0, 0, "cm"),
            text = element_text(size=9, family="Times"),
            panel.grid.minor = element_blank(),  # switch off minor gridlines
            panel.grid.major = element_blank(),
            legend.position.inside = c(0.78, 0.15), legend.justification = c(0, 1),
            axis.text = element_text(size=9, family="Times", color="black"),
            axis.title.y = element_text( angle=90),
            strip.text = element_text(size=9))
    #show(p_rec_plt)
    ggsave(p_rec_plt,
           file= paste0("figures/",cur_model, cur_ntrials, ".jpg"),width=15, height=15, units = "cm", dpi=600)
    ggsave(p_rec_plt,
           file= paste0("figures/",cur_model, cur_ntrials, ".eps"),
           width=17, height=19, units = "cm", dpi=600,
           device=cairo_ps)
  }
}


fitting_times <- pars %>% filter(class=="est") %>%
  select(model, fittingmins, ntrials)
ggplot(fitting_times, aes(x=as.factor(ntrials*10), y=fittingmins/60))+
  geom_boxplot()+
  facet_grid2(.~model, independent = "y", scales = "free_y")+
  theme_bw()+
  labs(y="Fitting time [hours]", x="Simulated trials")+
  theme(plot.margin = margin(0, 0.5, 0, 0, "cm"),
        text = element_text(size=12, family="Times"),
        panel.grid.major = element_blank(),
        legend.position.inside = c(0.78, 0.15), legend.justification = c(0, 1),
        axis.text = element_text(size=12, family="Times", color="black"),
        axis.title.y = element_text( angle=90),
        strip.text = element_text(size=12))
ggsave(file = "figures/recovery_fitting_times.eps", width = 18, height=7, 
       units = "cm", dpi=600, device = cairo_ps)
# ## Only for manuscript
# ggsave(file = "../../Draft/figures/recovery/recovery_fitting_times.eps",
#        width = 19, height=6.3,
#        units = "cm", dpi=600, device = cairo_ps)

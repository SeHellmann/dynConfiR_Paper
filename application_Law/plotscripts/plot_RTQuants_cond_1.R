windowsFonts(Times=windowsFont("Times"))
dir.create("figures", showWarnings = FALSE)

### Data response time quantiles                      ----  
# Reaction Time Quantiles of the Data grouped by condition and accuracy  
Data_RTQuants_corr_cond <- Data %>%
  group_by(condition, correct) %>%
  reframe(p=c(.1,.3,.5,.7,.9), q = quantile(rt, probs = c(.1,.3,.5,.7,.9))) 

### Prediction response time quantiles                ----
model_levels <-  c("PCRMt", "IRMt", "dynaViTE", "dynWEV", "2DSD")
Preds_RTQuants_corr_cond <-   Preds_RTdens_corr_cond_rating %>%
  group_by(model, rt, correct, condition) %>%
  reframe(dens = sum(dens)) %>%
  PDFtoQuantiles() %>%
  left_join(summarise(group_by(Preds_RatingDist_corr_cond,
                               model, condition, correct),
                      p_correct=sum(p), .groups="drop"),
            by = c("model","correct", "condition"))  %>%
  mutate(model=factor(model, levels=model_levels, ordered = TRUE))




## Figure 7: RTQuantiles accross correct X condition          ----
Data_RTQuants_corr_cond_plot <- Data_RTQuants_corr_cond %>%
  mutate(X = factor(paste(condition, correct, sep="."),
                    levels = c(paste(c(5:1, 1:5), rep(0:1, each=5), sep=".")),
                    labels = c(paste(c(5:1, 1:5), rep(0:1, each=5), sep="."))))%>% 
  filter(p %in% c(0.1, 0.5, 0.9))
Preds_RTQuants_corr_cond_plot <- Preds_RTQuants_corr_cond %>%
  mutate(X = factor(paste(condition, correct, sep="."),
                    levels = c(paste(c(5:1, 1:5), rep(0:1, each=5), sep=".")),
                    labels = c(paste(c(5:1, 1:5), rep(0:1, each=5), sep="."))),
         )%>% #model=str_replace(model, "\n" , " ") 
  filter(p %in% c(0.1, 0.5, 0.9))
ggplot()+
  geom_line(data=mutate(Preds_RTQuants_corr_cond_plot, correct=factor(correct, labels=c("Wrong", "Correct")),
                        condition = as.factor(condition)),
            aes(x=condition, y=log(q), group=as.factor(p),color=correct), linewidth=0.7)+
  geom_point(data=mutate(Data_RTQuants_corr_cond_plot, correct=factor(correct, labels=c("Wrong", "Correct")),
                         condition = as.factor(condition)),
             aes(x=condition, y=log(q), shape=correct),
             size=1.2, fill="white")+
  scale_color_manual(values= two_colors_correct, breaks=c("Correct", "Wrong"),
                     name = "Predicted",
                     labels=c("Correct", "Wrong")) +
  scale_shape_manual(values=c(21,17),breaks=c("Correct", "Wrong"),
                     name = "Observed",
                     labels=c("Correct", "Wrong"))  +
  scale_x_discrete(name="Motion coherence level", breaks=1:5)+
  scale_y_continuous(breaks = log(c(0.14, 0.3, 0.65, 1.4)),
                     # labels = paste("log(", c(1.5, 2.5, 4, 7), ")", sep=""), 
                     labels = format(c(0.14, 0.3, 0.65, 1.4), nsmall=2), 
                     name="Reaction time quantiles [s] (log scaled)")+
  facet_nested(model ~correct)+ #correct~model
  theme_bw() +
  theme(plot.margin = margin(0, 0, 0, 0, "cm"),
        text = element_text(size=9, family="Times"),
        legend.text = element_text(size=9),
        axis.text = element_text(size=9, family="Times", color="black"),
        axis.text.x = element_text(size=7, family="Times", color="black"),
        axis.title.y = element_text( angle=90),
        panel.grid.minor = element_blank(),  # switch off minor gridlines
        panel.grid.major = element_blank(),
        strip.text = element_text(size=9, family = "Times", angle=0),
        strip.text.y.right = element_text(angle = -90),
        legend.box = "horizontal",
        legend.position = "bottom",# legend.position = c(.70, .005), legend.justification = c(0,0),
        legend.direction = "horizontal", legend.spacing.y = unit(0, "lines"),
        legend.margin =margin(0,0,0,1, "cm"), legend.box.spacing = unit(0.2,"lines"),
        #legend.title = element_blank(),
        legend.key = element_blank(),
        legend.key.width=unit(1.5,"line"),
        panel.spacing=unit(0, "lines"))

ggsave("figures/RTQuantsCond_1.eps",
       width=17.62, height=15, dpi=1200, units="cm", device=cairo_ps)
ggsave("../../Draft/figures/example/RTQuantsCond_1.eps",
       width=17.62, height=15, dpi=1200, units="cm", device=cairo_ps)
ggsave("figures/RTQuantsCond_1.tiff",
       width = 15, height=13, units="cm",dpi=600)   # Filling a whole power point slide
# ggsave("../../Draft/figures/results/RTQuantsConf1.eps",
#        width=17.62, height=22, dpi=1200, units="cm", device=cairo_ps)
# ggsave("../../Draft/figures/results/RTQuantsConf1.eps",
#        width=18.288, height=15.24, dpi=1200, units="cm", device=cairo_ps)







             














# 
# 
# 
# ### Data response time quantiles                      ----  
# # Reaction Time Quantiles of the Data grouped by condition and accuracy  
# Data_RTQuants_corr_cond <- Data %>%
#   group_by(condition, correct, participant) %>%
#   reframe(p=c(.1,.3,.5,.7,.9), q = quantile(rt, probs = c(.1,.3,.5,.7,.9))) %>%
#   ungroup() %>%
#   Rmisc::summarySEwithin(measurevar = "q", idvar = "participant",
#                          withinvars = c("condition", "correct", "p"))
# 
# ### Prediction response time quantiles                ----
# Preds_RTQuants_corr_cond <-   prediction_RTConfDist %>%
#   group_by(participant, model, rt, condition, correct) %>%
#   reframe(dens=sum(dens)/2) %>%
#   group_by(model, participant, correct, condition) %>%
#   do(PDFtoQuantiles(.)) %>%
#   ungroup() %>%
#   group_by(model, condition, correct, p) %>%
#   reframe(q = mean(q))
# 
# 
# 
# 
# 
# ## Figure 7: RTQuantiles accross correct X condition          ----
# Data_RTQuants_corr_cond_plot <- Data_RTQuants_corr_cond %>%
#   filter(p %in% c(0.1, 0.5, 0.9))
# Preds_RTQuants_corr_cond_plot <- Preds_RTQuants_corr_cond %>%
#   filter(p %in% c(0.1, 0.5, 0.9)) %>%
#   mutate(correct=as.factor(correct), condition=as.factor(condition))
# ggplot(data=mutate(Data_RTQuants_corr_cond_plot, correct=factor(correct, labels=c("Wrong", "Correct")),
#                    condition = as.factor(condition)))+
#   geom_line(data=mutate(Preds_RTQuants_corr_cond_plot, correct=factor(correct, labels=c("Wrong", "Correct")),
#                         condition = as.factor(condition)),
#             aes(x=condition, y=log(q), group=as.factor(p),color=correct), linewidth=0.7)+
#   geom_point(aes(x=condition, y=log(q), shape=correct),
#              size=1.2, fill="white")+
#   geom_errorbar(aes(x=condition, ymin=log(q-se), ymax=log(q+se)), width=0.1)+
#   scale_color_manual(values= two_colors_correct, breaks=c("Correct", "Wrong"),
#                      name = "Predicted",
#                      labels=c("Correct", "Wrong")) +
#   scale_shape_manual(values=c(21,17),breaks=c("Correct", "Wrong"),
#                      name = "Observed",
#                      labels=c("Correct", "Wrong"))  +
#   scale_x_discrete(name="Motion coherence level", breaks=1:5)+
#   scale_y_continuous(breaks = log(c(0.14, 0.3, 0.65, 1.4)),
#                      # labels = paste("log(", c(1.5, 2.5, 4, 7), ")", sep=""), 
#                      labels = format(c(0.14, 0.3, 0.65, 1.4), nsmall=2), 
#                      name="Reaction time quantiles [s] (log scaled)")+
#   facet_nested(model ~correct)+ #correct~model
#   theme_bw() +
#   theme(plot.margin = margin(0, 0, 0, 0, "cm"),
#         text = element_text(size=9, family="Times"),
#         legend.text = element_text(size=9),
#         axis.text = element_text(size=9, family="Times", color="black"),
#         axis.text.x = element_text(size=7, family="Times", color="black"),
#         axis.title.y = element_text( angle=90),
#         panel.grid.minor = element_blank(),  # switch off minor gridlines
#         panel.grid.major = element_blank(),
#         strip.text = element_text(size=9, family = "Times", angle=0),
#         strip.text.y.right = element_text(angle = -90),
#         legend.box = "horizontal",
#         legend.position = "bottom",# legend.position = c(.70, .005), legend.justification = c(0,0),
#         legend.direction = "horizontal", legend.spacing.y = unit(0, "lines"),
#         legend.margin =margin(0,0,0,1, "cm"), legend.box.spacing = unit(0.2,"lines"),
#         #legend.title = element_blank(),
#         legend.key = element_blank(),
#         legend.key.width=unit(1.5,"line"),
#         panel.spacing=unit(0, "lines"))
# 
# ggsave("figures/RTQuantsCond.eps",
#        width=17.62, height=15, dpi=1200, units="cm", device=cairo_ps)
# ggsave("figures/RTQuantsCond.tiff",
#        width = 15, height=13, units="cm",dpi=600)   # Filling a whole power point slide
# # ggsave("../../Draft/figures/results/RTQuantsConf1.eps",
# #        width=17.62, height=22, dpi=1200, units="cm", device=cairo_ps)
# # ggsave("../../Draft/figures/results/RTQuantsConf1.eps",
# #        width=18.288, height=15.24, dpi=1200, units="cm", device=cairo_ps)
# 
# 
# 
# 
# 
# 

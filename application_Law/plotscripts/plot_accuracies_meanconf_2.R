windowsFonts(Times=windowsFont("Times"))
dir.create("figures", showWarnings = FALSE)

###########    Plot accuracies    #######
model_levels_2 <- c("dynaViTE\n(fixed)","dynaViTE","dynWEV\n(fixed)","dynWEV", "2DSD\n(fixed)", "2DSD")
model_colors_2<-c(  "#F0E442",      "#D55E00" ,     "#56B4E9",    "#0072B2","#000000",  "#CC79A7")
# model_levels_2 <- c("PCRMt",   "IRMt", "dynaViTE\n(fixed)","dynaViTE","dynWEV\n(fixed)","dynWEV", "2DSD\n(fixed)", "2DSD")
# model_colors_2<-c(   "#E69F00", "#009E73","#F0E442",      "#D55E00" ,     "#56B4E9",    "#0072B2","#000000",  "#CC79A7")

Preds_Acc <- rbind(prediction_ConfDist_fixed, prediction_ConfDist) %>% 
  group_by(participant, model, condition) %>%
  mutate(model= str_replace(model, " ", "\n")) %>%
  summarise(Acc = sum(p*correct)/(2))%>% 
  group_by(model, condition) %>%
  summarise(Acc = mean(Acc)) %>%
  filter(!grepl("RM",model)) %>%
  mutate(model=factor(model, levels=model_levels_2, ordered = TRUE)) 
Data_Acc <- Data %>% group_by(participant, condition) %>%
  summarise(Acc = mean(correct), .groups="drop") %>%
  reframe(Rmisc::summarySEwithin(pick(everything()),measurevar = "Acc", idvar = "participant", withinvars = c("condition"))) %>%
  rename(sewithin = se) 




## Figure: Plot of Fitted Accuracy                    ----
p_Acc <- ggplot(Data_Acc,
                aes(x=condition, y=Acc)) +
  geom_line(data=Preds_Acc, aes(x = condition, linetype="Predicted", group=model), linewidth=1)+
  #geom_line(linetype="dashed", alpha=0.5,position = pd)+
  geom_errorbar(aes(ymin=Acc-sewithin, ymax=Acc+sewithin), colour="black", width=.25, linewidth=0.4) +
  geom_point(fill="white", aes(shape="Observed"))+
#  facet_grid(cols=vars(model), drop=TRUE)+ #, dir="v"
  facet_wrap(.~model, nrow=1)+
  ylab("Mean accuracy")+
  scale_x_discrete(name="")+
  guides(x="none")+
  scale_linetype_manual(name="", values=1) +
  scale_shape_manual(values=c(21),name = "")  +
  scale_y_continuous(breaks=c(0.6, 0.7, 0.8), name="Mean Accuracy")+
  theme_bw() +
  guides(shape=guide_legend(order=3), color=guide_legend(order=3))+
  theme(plot.margin = margin(0, 0, 0, 0, "cm"),
        #axis.text.x = element_text(angle=90,hjust = 0.5, vjust = 0.5),
        text = element_text(size=9, family="Times"),
        axis.text = element_text(size=9, family="Times", color="black"),
        axis.title.y = element_text( angle=90),
        panel.grid.minor = element_blank(),  # switch off minor gridlines
        panel.grid.major = element_blank(),
        strip.text = element_text(size=9),
        legend.position = "right",
        # legend.position=c(.80, .005), 
        legend.margin = margin(0, 0.1, 0, 0, "cm"),#legend.direction = "horizontal",
        legend.box.margin = margin(-0.1, 0, 0, 0, "cm"),
        legend.key.width=unit(1,"line"),
        legend.key.height = unit(0.2, "line"), legend.spacing.y = unit(0.2, "line"),
        panel.spacing=unit(0, "lines"))
p_Acc







Data_MRating_corr_cond_part <- Data_RatingDist_part %>% group_by(participant, condition, correct) %>%
  summarise(meanRating = sum(p*rating)/sum(p), .groups = "drop") %>%
  mutate(condition=as.factor(condition), correct=as.factor(correct))
Data_MRating_corr_cond <- Data_MRating_corr_cond_part %>%
  Rmisc::summarySEwithin(measurevar = "meanRating", withinvars = c("condition", "correct"), 
                         idvar = "participant")
Preds_MRating_corr_cond_part <- rbind(prediction_ConfDist, prediction_ConfDist_fixed) %>% 
  group_by(model, participant, condition, correct) %>%
  summarise(meanRating = sum(p*rating)/(sum(p)), .groups = "drop") %>%
  mutate(condition=as.factor(condition), correct=as.factor(correct))
Preds_MRating_corr_cond <- Preds_MRating_corr_cond_part %>%
  mutate(model= str_replace(model, " ", "\n")) %>%
  group_by(model, condition, correct) %>%
  summarise(meanRating=mean(meanRating)) %>%
  filter(!grepl("RM",model)) %>%
  mutate(model=factor(model, levels=model_levels_2, ordered=TRUE))

two_colors_correct <- c("#1b9e77", "#fc8d62")

## Figure 6: Plot of Mean Ratings                          ----
Preds_MRating_corr_cond_plot <- mutate(Preds_MRating_corr_cond,
                                       confidence ="Confidence measure")
Data_MRating_corr_cond_plot <- merge(Data_MRating_corr_cond, 
                                     expand.grid(confidence ="Confidence measure"))
p_MRating <- ggplot(Data_MRating_corr_cond_plot,
                    aes(x=condition, y=meanRating, group = as.factor(correct), shape=as.factor(correct))) +
  # geom_ribbon(data=Preds_MRating_corr_cond, 
  #             aes(ymin=meanRating-sewithin, ymax=meanRating+sewithin, fill=as.factor(correct)), alpha=0.5)+ #
  geom_line(data=Preds_MRating_corr_cond_plot, aes(color=as.factor(correct)), linewidth=0.8)+
  #geom_line(linetype="dashed", alpha=0.5,position = pd)+
  geom_errorbar(aes(ymin=meanRating-se, ymax=meanRating+se), colour="black", width=.25,linewidth=0.6) +
  geom_point(fill="white", size=1.8)+
  # facet_nested(rows=vars(model), drop=TRUE)+ #, dir="v"
  facet_wrap(.~model, nrow=1)+ #, dir="v"
  scale_x_discrete(name="Motion coherence level")+
  scale_color_manual(values= two_colors_correct, breaks=c(1,0),
                     name = "Predicted", labels=c("Correct", "Wrong")) +
  scale_fill_manual(values= two_colors_correct, breaks=c(1,0),
                    name = "Predicted", labels=c("Correct", "Wrong")) +
  scale_shape_manual(values=c(21,17),breaks=c(1,0),
                     name = "Observed", labels=c("Correct", "Wrong"))  +
  scale_y_continuous(breaks = c(2.0, 2.3, 2.6),name="Mean confidence rating")+
  theme_bw() +
  guides(shape=guide_legend(order=3), color=guide_legend(order=3), fill="none")+
  theme(plot.margin = margin(-0.1, 0, 0, 0, "cm"),
        text = element_text(size=9, family="Times"),
        axis.text = element_text(size=9, family="Times", color="black"),
        axis.title.y = element_text( angle=90),
        panel.grid.minor = element_blank(),  # switch off minor gridlines
        panel.grid.major = element_blank(),
        strip.text = element_blank(), 
        strip.background = element_blank(),
        legend.position = "right", 
        legend.margin = margin(0,0.2,0,0.1,"cm"),
        legend.direction = "vertical", legend.justification = c(0,0),
        legend.box="vertical",
        legend.key = element_blank(), legend.spacing = unit(0,"line"),
        legend.key.width=unit(1,"line"), 
        legend.title = element_text(margin=margin(t=1, unit="line")),
        legend.key.height = unit(0, "line"), legend.spacing.y = unit(0, "line"),
        legend.text=element_text(margin=margin(b=0.2, unit="cm")),
        panel.spacing=unit(0, "lines"))
p_MRating
ggarrange(p_Acc, p_MRating, ncol=1, heights = c(0.48, 0.52))
dir.create("figures", showWarnings = FALSE)
ggsave("figures/accuracymeanRating_2.tiff",
       width = 17.62, height=12, units="cm",dpi=600)
ggsave("figures/accuracymeanRating_2.eps",
       width = 15, height=7, units="cm",dpi=1200, device = cairo_ps)
# ggsave("../../Draft/figures/example/accuracymeanRating_2.eps",
#        width = 15, height=7, units="cm",dpi=1200, device = cairo_ps)



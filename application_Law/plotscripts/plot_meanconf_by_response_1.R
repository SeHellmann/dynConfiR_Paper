windowsFonts(Times=windowsFont("Times"))
dir.create("figures", showWarnings = FALSE)


Data_MRating_corr_cond_resp_part <- Data_RatingDist_part %>% 
  group_by(participant, condition, correct, response) %>%
  summarise(meanRating = sum(p*rating)/sum(p), .groups = "drop") %>%
  mutate(condition=as.factor(condition), correct=as.factor(correct)) %>%
  filter(!is.na(meanRating))
Data_MRating_corr_resp_cond <- Data_MRating_corr_cond_resp_part %>%
  mutate(response=factor(response, levels=c(1,2), labels=c("Leftward motion", "Rightward motion"))) %>%
  Rmisc::summarySEwithin(measurevar = "meanRating", withinvars = c("condition", "correct", "response"), 
                         idvar = "participant")
Preds_MRating_corr_cond_resp_part <- prediction_ConfDist %>% group_by(model, participant, condition, correct, response) %>%
  summarise(meanRating = sum(p*rating)/(sum(p)), .groups = "drop") %>%
  mutate(response=ifelse(!grepl("RM", model), response*0.5+1.5, response)) %>%
  mutate(condition=as.factor(condition), correct=as.factor(correct),
         response=factor(response, levels=c(1,2), labels=c("Leftward motion", "Rightward motion")))
Preds_MRating_corr_resp_cond <- Preds_MRating_corr_cond_resp_part %>%
  group_by(model, condition, correct, response) %>%
  summarise(meanRating=mean(meanRating)) %>%
  mutate(model=factor(model, levels=model_levels, ordered=TRUE))

two_colors_correct <- c("#1b9e77", "#fc8d62")

## Figure 6: Plot of Mean Ratings                          ----
Preds_MRating_corr_cond_resp_plot <- mutate(Preds_MRating_corr_resp_cond,
                                       confidence ="Confidence measure")
Data_MRating_corr_cond_resp_plot <- merge(Data_MRating_corr_resp_cond, 
                                     expand.grid(confidence ="Confidence measure"))
p_MRating_resp <- ggplot(Data_MRating_corr_cond_resp_plot,
                    aes(x=condition, y=meanRating, group = as.factor(correct), shape=as.factor(correct))) +
  # geom_ribbon(data=Preds_MRating_corr_cond, 
  #             aes(ymin=meanRating-sewithin, ymax=meanRating+sewithin, fill=as.factor(correct)), alpha=0.5)+ #
  geom_line(data=Preds_MRating_corr_cond_resp_plot, aes(color=as.factor(correct)), linewidth=0.8)+
  #geom_line(linetype="dashed", alpha=0.5,position = pd)+
  geom_errorbar(aes(ymin=meanRating-se, ymax=meanRating+se), colour="black", width=.25,linewidth=0.6) +
  geom_point(fill="white", size=1.8)+
  # facet_nested(rows=vars(model), drop=TRUE)+ #, dir="v"
  facet_grid(cols=vars(response), rows=vars(model))+ #, dir="v"
  scale_x_discrete(name="Motion coherence level")+
  scale_color_manual(values= two_colors_correct, breaks=c(1,0),
                     name = "Predicted", labels=c("Correct", "Wrong")) +
  scale_fill_manual(values= two_colors_correct, breaks=c(1,0),
                    name = "Predicted", labels=c("Correct", "Wrong")) +
  scale_shape_manual(values=c(21,17),breaks=c(1,0),
                     name = "Observed", labels=c("Correct", "Wrong"))  +
  scale_y_continuous(breaks = c(2.0, 2.4, 2.8),name="Mean confidence rating")+
  theme_bw() +
  guides(shape=guide_legend(order=3), color=guide_legend(order=3), fill="none")+
  theme(plot.margin = margin(-0.1, 0, 0, 0, "cm"),
        text = element_text(size=9, family="Times"),
        axis.text = element_text(size=9, family="Times", color="black"),
        axis.title.y = element_text( angle=90),
        panel.grid.minor = element_blank(),  # switch off minor gridlines
        panel.grid.major = element_blank(),
        # strip.text = element_blank(), 
        # strip.background = element_blank(),
        legend.position = "bottom", 
        legend.margin = margin(0,0.2,0,0.1,"cm"),
        # legend.direction = "vertical", legend.justification = c(0,0),
        legend.box="vertical",
        legend.key = element_blank(), legend.spacing = unit(0,"line"),
        legend.key.width=unit(1,"line"), 
        #legend.title = element_text(margin=margin(t=1, unit="line")),
        legend.key.height = unit(0, "line"), legend.spacing.y = unit(0, "line"),
        legend.text=element_text(margin=margin(b=0.2, unit="cm")),
        panel.spacing=unit(0, "lines"))
p_MRating_resp
dir.create("figures", showWarnings = FALSE)
ggsave("figures/meanRating_response_1.tiff",
       width = 17.62, height=15, units="cm",dpi=600)
ggsave("figures/meanRating_response_1.eps",
       width = 12, height=15, units="cm",dpi=1200, device = cairo_ps)
ggsave("../../Draft/figures/example/meanRating_response_1.eps",
       width = 12, height=15, units="cm",dpi=1200, device = cairo_ps)



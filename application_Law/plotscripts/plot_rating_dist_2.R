Data <- Data %>% group_by(participant, condition, stimulus) %>%
  mutate(nrows=n())
Data_RatingDist_part <- Data %>% 
  group_by(participant, condition, rating, correct, stimulus) %>% 
  summarise(p = n()/(mean(nrows)),.groups = "drop") %>%
  full_join(y = expand.grid(rating=1:4, correct=0:1, stimulus=c(1,2),
                            condition=1:5, participant=unique(Data$participant)),
            by = c("participant", "condition", "rating", "correct", "stimulus")) %>%
  mutate(p = ifelse(is.na(p), 0, p)) 

### Sanity Checks:  15 (participants)*5(conditions):
#sum(Data_RatingDist_part$p)
#table((Data_RatingDist_part %>% group_by(condition, participant) %>% summarise(p = sum(p)))$p)

# For the plots we won't differentiate between 
# stimulus directions and participants
Data_RatingDist_corr_cond <- Data_RatingDist_part %>% 
  group_by(correct, condition, rating, stimulus) %>% 
  summarise(p = mean(p), .groups = "drop")

### d) Prediction confidence rating distribution         ----
# model_levels_2 <- c("PCRMt",   "IRMt", "dynaViTE\n(fixed)","dynaViTE","dynWEV\n(fixed)","dynWEV", "2DSD\n(fixed)", "2DSD")
# model_colors_2<-c(   "#E69F00", "#009E73","#F0E442",      "#D55E00" ,     "#56B4E9",    "#0072B2","#000000",  "#CC79A7")
# shapes = c(           22,          24,        23,             23,             21,          21,        25,       25)

model_levels_2 <- c("dynaViTE\n(fixed)","dynaViTE","dynWEV\n(fixed)","dynWEV", "2DSD\n(fixed)", "2DSD")
model_colors_2<-c(  "#F0E442",      "#D55E00" ,     "#56B4E9",    "#0072B2","#000000",  "#CC79A7")
shapes = c(         23,             23,             21,          21,        25,       25)
Preds_RatingDist_corr_cond <- rbind(prediction_ConfDist_fixed, prediction_ConfDist) %>% 
  mutate(model= str_replace(model, " ", "\n")) %>%
  group_by(model, rating, correct, condition, stimulus) %>%
  mutate(stimulus=ifelse(!grepl("RM", model), stimulus*0.5+1.5, stimulus)) %>%
  filter(!grepl("RM", model)) %>%   ##### RM raus
  summarise(p = mean(p), .groups="drop") %>%
  mutate(model=factor(model, levels=model_levels_2, ordered=TRUE))
## More Sanity checks:
# Data_RatingDist_corr_cond %>% group_by(condition) %>% summarise(p=sum(p))
# Preds_RatingDist_corr_cond %>% group_by(condition, model) %>% summarise(p=sum(p))

## Suppl Fig 7-9: Fitted response distributions            ----
# for the colors: https://davidmathlogic.com/colorblind/#%23D81B60-%231E88E5-%23FFC107-%23004D40 
p_ratingdist <- ggplot(data=Data_RatingDist_corr_cond, aes(x=rating, y=p))+
  geom_bar(stat = "identity", show.legend = FALSE, fill="white", col="black")+
  geom_point(data=Preds_RatingDist_corr_cond,
             aes(shape=model, fill=model), 
             position=position_dodge(1))+
  facet_nested(cols=c(vars(stimulus), vars(correct)), rows=c(vars("Coherence level"),vars(condition)),
               labeller = labeller(correct=c("0"="Wrong", "1"="Correct"),
                                   stimulus=c("1"="Leftwards motion", "2"="Rightwards motion")),
               strip=strip_nested(by_layer_y = TRUE, text_y = list(element_text(angle = -90), NULL)))+
  scale_x_continuous(name = "Confidence level", breaks = 1:4)+
  scale_y_continuous(name = "Probability", breaks = c(0, 0.1,0.2, 0.3))+
  scale_shape_manual(name = "Model prediction",values=shapes) +
  scale_fill_manual(name = "Model prediction",values=model_colors_2) +
  theme_bw() +
  theme(plot.margin = margin(0, 0.2, 0, 0, "cm"),
        text = element_text(size=9, family="Times"),
        axis.text = element_text(size=9, family="Times", color="black"),
        #axis.text.x = element_text(angle=90),
        axis.title.y = element_text( angle=90),
        panel.grid.minor = element_blank(),  # switch off minor gridlines
        panel.grid.major = element_blank(),
        legend.key = element_blank(),
        legend.position = "bottom", legend.text = element_text(size=9, family = "Times"),
        strip.text.y = element_text(angle = 0, size=9, family = "Times"),
        panel.spacing=unit(0, "lines"))+
  coord_cartesian(xlim = c(0.5,4.5),ylim=c(0, 0.35),  clip = 'off')
p_ratingdist <- p_ratingdist+
  theme(plot.margin = margin(0, 0, 0, 0, "cm"),
        legend.box = "vertical", legend.spacing.y = unit(-.2, "cm"))
p_ratingdist
ggsave(paste0("figures/RespDist_2.tiff"), plot = p_ratingdist,
     width=17.3, height=18, units="cm",dpi=1200)
ggsave("figures/RespDist_2.eps",
     width=17.3, height=18, units="cm",dpi=1200, device = cairo_ps)
ggsave("../../Draft/figures/example/RespDist_2.eps",
       width=19, height=18, units="cm",dpi=1200, device = cairo_ps)





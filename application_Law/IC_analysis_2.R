require(ggpubr)
require(scales)
windowsFonts(Times=windowsFont("Times"))
dir.create("figures", showWarnings = FALSE)

temp_allfits <- allfits %>% 
  mutate(model= str_replace(model, " ", "\n")) %>% 
  select(model, participant, BIC, AIC, AICc) %>%#filter(!grepl("DDM", model)) %>%
  pivot_longer(cols=c("BIC", "AIC", "AICc"), names_to="criterion") %>%
  mutate(model = as.factor(model), criterion=as.factor(criterion))

plotICs <- temp_allfits %>%
  group_by(model, criterion) %>%
  summarise(Mean=mean(value), SE=sd(value)/sqrt(n()), .groups = "drop")
plotICs <- temp_allfits %>%
  group_by(criterion) %>%
  reframe(Rmisc::summarySEwithin(pick(everything()), measurevar="value", 
                                 idvar="participant", withinvars = "model")) %>%
  right_join(plotICs)

plotICs <- plotICs %>% 
  mutate(model= factor(model, 
                       levels=plotICs[plotICs$criterion=="BIC",][["model"]][order(plotICs[plotICs$criterion=="BIC",][["Mean"]])], ordered = TRUE))

ggplot(temp_allfits, aes(x=model, y=value))+
  geom_violin()+geom_line(aes(group=participant))+
  facet_wrap(.~criterion)

temp_allfits2 <- temp_allfits %>% group_by(criterion, participant) %>%
  mutate(value = value - mean(value), 
         class="Mean-centered") %>%
  rbind(mutate(temp_allfits, class="Absolute values"))  %>%
  ungroup()
ggplot(subset(temp_allfits2,criterion =="BIC"), aes(x=model, y=value))+
  geom_violin()+geom_line(aes(group=participant))+
  facet_grid(class~., scales = "free_y")

ggplot(subset(temp_allfits2, criterion =="BIC"), aes(x=model, y=value))+
  geom_violin()+geom_line(aes(group=participant))+
  facet_grid(class~., scales = "free_y")
distinct(select(temp_allfits, model, participant, criterion))
temp2 <- temp_allfits %>% 
  filter(!grepl("dynaViTE", model) ) %>%
  left_join(select(filter(temp_allfits, grepl("dynaViTE", model)&grepl("fixed", model)), 
                   participant, criterion, value_dynaViTE=value)) %>%
  mutate(value = value-value_dynaViTE)

ggplot(subset(temp2,criterion=="BIC"), aes(x=model, y=value))+
  #geom_violin()+
  geom_line(aes(group=participant))+geom_point()+
  ylab("Difference in IC between model and dynaViTE (fixed)")
ggplot(subset(temp2,criterion=="BIC"), aes(x=model, y=sign(value)*log(abs(value), base=3)))+
  #geom_violin()+
  geom_line(aes(group=participant))+geom_point()+
  scale_y_continuous(breaks=-5:5)+
  ylab("Difference in IC between model and dynaViTE (fixed)")

### 
model_levels_2 <- levels(plotICs$model)
model_levels_2
# Colorblind-friendly color palette:
model_levels <-  c(          "PCRMt", "IRMt",              "dynaViTE", "dynWEV", "2DSD")
model_colors <- c(          "#E69F00","#009E73",            "#D55E00", "#0072B2", "#CC79A7")
cbbPalette <- c("#F0E442", "#E69F00", "#009E73", "#56B4E9", "#D55E00" ,"#0072B2",  "#CC79A7","#000000")

# To use for fills, add
p_BICavg <- ggplot(subset(plotICs,criterion=="BIC"), 
                   aes(group=1, x=model, y=-Mean, ymin=-Mean-se, ymax=-Mean+se))+
  #geom_line(linewidth=2)+
  geom_point(aes(color=as.character(model)), shape=16, size=4)+
  scale_color_manual(name="Model", 
                     breaks = levels(plotICs$model),
                     values=cbbPalette, guide="none")+
  geom_errorbar(aes(color=as.character(model)),width=0.2, linewidth=1)+
  scale_y_continuous(name = "Mean negative\nBIC (thousands)",
                     breaks = c(-3, -4, -5, -6, -7, -8,-10)*1000,
                     labels=label_number(suffix = "", scale = 1e-3))+
  scale_x_discrete(name="Model")+
  #facet_grid(cols=vars(criterion))+
  theme_bw()+
  theme(plot.margin = margin(0, 0, 0, 0.1, "cm"),
        text = element_text(size=12, family="Times"),
        panel.grid.minor = element_blank(),  # switch off minor gridlines
        panel.grid.major.x = element_blank(),
        axis.text = element_text(size=12, family="Times", color="black"),
        axis.text.x = element_text(angle=45, hjust=1),
        axis.title.y = element_text( angle=90))
p_BICavg

# ggsave("figures/ICs_2.eps", device=cairo_ps, width=8, height=3.5, dpi=600)
# ggsave("../../Draft/figures/example/ICs_2.eps", device=cairo_ps, 
#        width=8, height=3.5, dpi=600)
# 
# ggsave("figures/ICs_2.tiff", width=8, height=3.5, dpi=600)



BICweights <- allfits %>%
  mutate(model= str_replace(model, " ", "\n")) %>% 
  group_by(participant) %>%
  mutate(BICdiff = BIC-min(BIC),
         expBICdiff = exp(-0.5*(BIC-min(BIC))),
         wBIC = expBICdiff / sum(expBICdiff)) %>%
  #  wBIC = round((expBICdiff+1e-300) / sum(expBICdiff+1e-300),20)) %>%
  ungroup() 
## Rearrange plot order for participants
BICweights <- BICweights %>% 
  ungroup() %>% 
  select(model,participant, wBIC) %>% 
  pivot_wider(names_from = model, values_from = wBIC) %>%
  arrange(round(desc(`dynaViTE`)),desc(`dynaViTE\n(fixed)`), desc(`2DSD`), desc(IRMt), desc(PCRMt)) %>%
  mutate(plotorder = 1:n()) %>% 
  pivot_longer(2:9, names_to="model", values_to="wBIC")
# BICweights <- BICweights %>% arrange(model, participant)
# BICweights <- BICweights %>%
#   group_by(model) %>% mutate(plotorder = order(BICweights[BICweights$model=="dynaViTE\n(fixed)",]$wBIC,
#                                                decreasing = TRUE))
BICweights <- mutate(BICweights, model=factor(model, levels=model_levels_2))
p_BIC_part<- ggplot(BICweights, aes(x=plotorder, y=wBIC, fill=model))+
  geom_bar(stat = "identity")+
  scale_fill_manual(name="Model",
                    breaks = levels(plotICs$model)[c(5,1,7, 3, 8, 2, 4, 6)],
                    values=cbbPalette[c(5,1,7, 3, 8, 2, 4, 6)])+
  ylab("BIC weight")+
  scale_x_continuous(name="Participant (reordered)", breaks=c(1, 5, 10, 15))+
  theme_minimal()+
  theme(plot.margin = margin(0.1, 0, 0, 0.1, "cm"),
        text = element_text(size=12, family="Times"),
        panel.grid.minor = element_blank(),  # switch off minor gridlines
        panel.grid.major.x = element_blank(),
        axis.text = element_text(size=12, family="Times", color="black"),
        axis.title.y = element_text( angle=90),
        strip.text = element_blank(), strip.background = element_blank(),
        legend.key.height = unit(0.7, "cm"),
        legend.key.width = unit(0.7, "cm"))
p_BIC_part
# ggsave("figures/BICweights_2.eps", width=15, height=6, 
#        units = "cm", device = cairo_ps)
# ggsave("../../Draft/figures/example/BICweights_2.eps", width=15, height=6, 
#        units = "cm", device = cairo_ps)
# ggsave("figures/BICweights_2.tiff", width=8.2, height=4, 
#        units = "cm")


ggarrange(p_BICavg, p_BIC_part, nrow=2, heights = c(0.6, 0.4))
ggsave("figures/BIC_2.eps", width=15, height=13, 
       units = "cm", device = cairo_ps)
ggsave("../../Draft/figures/example/BIC_2.eps", width=15, height=13, 
       units = "cm", device = cairo_ps)
ggsave("figures/BIC_2.tiff", width=8.2, height=4, 
       units = "cm")



require(ggpubr)
require(scales)
windowsFonts(Times=windowsFont("Times"))
dir.create("figures", showWarnings = FALSE)

temp_parfits <- parfits %>% 
  select(model, participant, BIC, AIC, AICc) %>%#filter(!grepl("DDM", model)) %>%
  pivot_longer(cols=c("BIC", "AIC", "AICc"), names_to="criterion") %>%
  mutate(model = as.factor(model), criterion=as.factor(criterion))

plotICs <- temp_parfits %>%
  group_by(criterion) %>%
  reframe(Rmisc::summarySEwithin(pick(everything()), measurevar="value", 
                                 idvar="participant", withinvars = "model")) %>%
  rename(Mean=value)

plotICs <- plotICs %>% 
  mutate(model= factor(model, 
                       levels=plotICs[plotICs$criterion=="BIC",][["model"]][order(plotICs[plotICs$criterion=="BIC",][["Mean"]])], ordered = TRUE))

# ggplot(temp_parfits, aes(x=model, y=value))+
#   geom_violin()+geom_line(aes(group=participant))+
#   facet_wrap(.~criterion)


### 
# Colorblind-friendly color palette:
cbbPalette <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
# To use for fills, add
p_BICavg <- ggplot(subset(plotICs,criterion=="BIC"), 
            aes(group=1, x=model, y=-Mean, ymin=-Mean-se, ymax=-Mean+se))+
  #geom_line(linewidth=2)+
  geom_point(aes(color=as.character(model)), shape=16, size=4)+
  scale_color_manual(name="Model", 
                     breaks = levels(plotICs$model),
                     values=cbbPalette[c(1, 3, 6, 5, 7)], guide="none")+
  geom_errorbar(aes(color=as.character(model)),width=0.2, linewidth=1)+
  scale_y_continuous(name = "Mean negative\nBIC (thousands)",
                     breaks = c(-3200, -3500, -3800, -4100, -4400),
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

ggsave("figures/ICs_1.eps", device=cairo_ps, width=8, height=3.5, dpi=600)
ggsave("../../Draft/figures/example/ICs_1.eps", device=cairo_ps, 
       width=8, height=3.5, dpi=600)

ggsave("figures/ICs_1.tiff", width=8, height=3.5, dpi=600)



BICweights <- parfits %>%
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
  arrange(desc(`dynaViTE`), desc(`2DSD`), desc(IRMt), desc(PCRMt)) %>%
  mutate(plotorder = 1:n()) %>% 
  pivot_longer(2:6, names_to="model", values_to="wBIC")
# BICweights <- BICweights %>% arrange(model, participant)
# BICweights <- BICweights %>%
#   group_by(model) %>% mutate(plotorder = order(BICweights[BICweights$model=="dynaViTE\n(fixed)",]$wBIC,
#                                                decreasing = TRUE))
p_BIC_part<- ggplot(BICweights, aes(x=plotorder, y=wBIC, fill=model))+
  geom_bar(stat = "identity")+
  scale_fill_manual(name="Model", 
                     breaks = levels(plotICs$model),
                     values=cbbPalette[c(1, 3, 6, 5, 7)])+
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
ggsave("figures/BICweights_1.eps", width=15, height=6, 
       units = "cm", device = cairo_ps)
ggsave("../../Draft/figures/example/BICweights_1.eps", width=15, height=6, 
       units = "cm", device = cairo_ps)
ggsave("figures/BICweights_1.tiff", width=8.2, height=4, 
       units = "cm")


ggarrange(p_BICavg, p_BIC_part, nrow=2, heights = c(0.6, 0.4))
ggsave("figures/BIC_1.eps", width=15, height=10, 
       units = "cm", device = cairo_ps)
ggsave("../../Draft/figures/example/BIC_1.eps", width=15, height=10, 
       units = "cm", device = cairo_ps)
ggsave("figures/BIC_1.tiff", width=8.2, height=4, 
       units = "cm")




plotICs$model
temp_parfits2 <- temp_parfits %>% group_by(criterion, participant) %>%
  mutate(value = value - mean(value), 
         class="Mean-centered") %>%
  ungroup() %>%
  mutate(model=factor(model, levels=levels(plotICs$model), ordered=T))
ggplot(subset(temp_parfits2,criterion=="BIC"), aes(x=model, y=value, fill=model))+
  geom_violin()+geom_line(aes(group=participant))+
  scale_fill_manual(name="Model", 
                    breaks = levels(plotICs$model),
                    values=cbbPalette[c(1, 3, 6, 5, 7)], guide="none")+
  xlab("Model")+ylab("BIC (centered for each participant)")+
  theme_bw()+
  theme(plot.margin = margin(0, 0, 0, 0.1, "cm"),
        text = element_text(size=12, family="Times"),
        panel.grid.minor = element_blank(),  # switch off minor gridlines
        panel.grid.major.x = element_blank(),
        axis.text = element_text(size=12, family="Times", color="black"),
        axis.title.y = element_text( angle=90))
ggsave("figures/CenteredBIC_1.eps", width=15, height=10, 
       units = "cm", device = cairo_ps)
ggsave("../../Draft/figures/example/CenteredBIC_1.eps", width=15, height=10, 
       units = "cm", device = cairo_ps)
ggsave("figures/CenteredBIC_1.tiff", width=15, height=10, 
       units = "cm")


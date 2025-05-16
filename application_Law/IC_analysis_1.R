require(ggpubr)
require(scales)
windowsFonts(Times=windowsFont("Times"))
dir.create("figures", showWarnings = FALSE)
custom_theme <- theme(plot.margin = margin(0.1, 0, 0, 0.1, "cm"),
        text = element_text(size=12, family="Times"),
        panel.grid.minor = element_blank(),  # switch off minor gridlines
        panel.grid.major.x = element_blank(),
        axis.text = element_text(size=12, family="Times", color="black"),
        axis.title.y = element_text( angle=90),
        strip.text = element_blank(), strip.background = element_blank(),
        legend.key.height = unit(0.7, "cm"),
        legend.key.width = unit(0.7, "cm"))
# Colorblind-friendly color palette:
cbbPalette <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")


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
p_BICavg <- ggplot(subset(plotICs,criterion=="BIC"), 
            aes(group=1, x=model, y=-Mean, ymin=-Mean-se, ymax=-Mean+se))+
  #geom_line(linewidth=2)+
  geom_point(aes(color=as.character(model)), shape=16, size=4)+
  scale_color_manual(name="Model", 
                     breaks = levels(plotICs$model),
                     values=cbbPalette[c(1, 3, 6, 5, 7)], guide="none")+
  geom_errorbar(aes(color=as.character(model)),width=0.2, linewidth=1)+
  scale_y_continuous(name = "Mean negative\nBIC (thousands)",
                     breaks = c(-3500, -4000, -4500, -5000),
                     labels=label_number(suffix = "", scale = 1e-3))+
  scale_x_discrete(name="Model")+
  #facet_grid(cols=vars(criterion))+
  theme_bw()+custom_theme+
  theme(axis.text.x = element_text(angle=45, hjust=1))
p_BICavg

### Group-level BMS
groupBMS <- group_BMS_fits(parfits)
groupBMS_df <- groupBMS$model_weights %>% 
  as.data.frame() %>%
  rownames_to_column("model") %>% 
  arrange(desc(pep)) %>% 
  mutate(pltorder = 1:n(), 
         model = factor(pltorder, labels=model))
# For dirichlet-rv we would scale Gamma-RV to a max of 1,
# but we are interested in the index of the max only, so this it not
# necessary
alpha <- groupBMS_df$alpha
n <- 4000
sim_dir <- expand.grid(N=1:n, model=unique(groupBMS_df$model))
sim_dir <- sim_dir %>% left_join(groupBMS_df[,c("model", "alpha", "fx_prob")]) %>%
  mutate(GammaRV = rgamma(1:n(), shape=alpha)) %>%
  group_by(N) %>% mutate(DirProb = GammaRV/sum(GammaRV)) %>% ungroup() %>%
  mutate(model=factor(model, levels=levels(plotICs$model)))

pep_plt <- ggplot(sim_dir)+
  geom_violin(aes(x=model, y=DirProb, fill=model))+
  xlab("Model")+ylab("Probability")+
  geom_errorbar(data=groupBMS_df, aes(x=model, ymin=pep, ymax=pep, linetype="Protected exceedance probability"), 
                width=0.35, linewidth=1.2)+#, shape="-", size="\U2014")+
  geom_errorbar(data=groupBMS_df, aes(x=model, ymin=fx_prob, ymax=fx_prob, linetype="Model probability (fixed effect model)"), 
                width=0.35, linewidth=1.2)+#, shape="-", size="\U2014")+
  scale_linetype_manual(name="", 
                        breaks=c("Protected exceedance probability", "Model probability (fixed effect model)"),
                        values=c(1,3))+
  scale_fill_manual(name="Model", 
                    breaks = levels(plotICs$model),
                    values=cbbPalette[c(1, 3, 6, 5, 7)], guide="none")+
  theme_bw()+custom_theme +
  theme(legend.position = "bottom",
        axis.text.x = element_text(angle=45, hjust=1),
        legend.key.width = unit(1, "cm"), legend.key.height = unit(0.1, "cm"))
pep_plt

BICweights <- subject_modelweights(parfits)
BICweights <- BICweights %>% 
  arrange(desc(`dynaViTE`), desc(`2DSD`), desc(IRMt), desc(PCRMt)) %>%
  mutate(plotorder = 1:n()) %>% 
  pivot_longer(cols=-c(participant, plotorder), names_to="model", values_to="wBIC") %>%
  mutate(model=factor(model, levels=levels(plotICs$model)))


p_BIC_part<- ggplot(BICweights, aes(x=plotorder, y=wBIC, fill=model))+
  geom_bar(stat = "identity", color="black", linewidth=0.6, width=1)+
  scale_fill_manual(name="Model", 
                     breaks = levels(plotICs$model),
                     values=cbbPalette[c(1, 3, 6, 5, 7)])+
  ylab("BIC weight")+
  scale_x_continuous(name="Participant (reordered)", breaks=c(1, 5, 10, 15), expand = c(0.02, 0.02))+
  theme_minimal()+custom_theme+
  theme(legend.position = "top")
p_BIC_part




p_BICavg2 <- p_BICavg + theme(axis.title.x = element_blank(), 
                              axis.text.x = element_blank())
p_BICavg2
ggarrange(p_BIC_part, p_BICavg2, pep_plt,nrow=3, heights = c(0.32, 0.22, 0.46), align = "v")
ggsave("figures/BIC_1.eps", width=15, height=19, 
       units = "cm", device = cairo_ps)
# ggsave("../../Draft/figures/example/BIC_1.eps", width=15, height=19,
#        units = "cm", device = cairo_ps)
ggsave("figures/BIC_1.tiff", width=15, height=19, 
       units = "cm")




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
# ggsave("../../Draft/figures/example/CenteredBIC_1.eps", width=15, height=10,
#        units = "cm", device = cairo_ps)
ggsave("figures/CenteredBIC_1.tiff", width=15, height=10, 
       units = "cm")


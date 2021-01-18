library(ggplot2)

noaa<-read.csv("P2S_CRED Data.csv")
biomass<-read.csv("Biomass_data_with_inhabitation.csv")
#______
#Overview

biomass_plot<-ggplot(biomass, aes(x = biomass$Calcifying, y = biomass$Microbial.Biomass.ug.l, shape= biomass$Inhabitation, fill = biomass$Inhabitation))+
  geom_point()+
  stat_smooth(method = "lm")+
  theme_bw()+
  xlab("Calcifying benthic cover (%)")+
  ylab("Microbial Biomass (\U00B5g/L)")+
  scale_shape_manual(values = c(21, 24))+
  scale_fill_manual(values = c("Grey", "Black"))+
  guides(fill=FALSE, shape = FALSE)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  ggsave("Calcifying x Biomass.pdf", width=8, height=8, units="cm")

biomass_plot

lm_biomass<- lm(biomass$Microbial.Biomass.ug.l ~ biomass$Calcifying)
summary(lm_biomass)

inhabited<-biomass[which(biomass$Inhabitation == "Inhabited"),]
lm_biomass_in<- lm(inhabited$Microbial.Biomass.ug.l ~ inhabited$Calcifying)
summary(lm_biomass_in)

uninhabited<-biomass[which(biomass$Inhabitation == "Uninhabited"),]
lm_biomass_un<- lm(uninhabited$Microbial.Biomass.ug.l ~ uninhabited$Calcifying)
summary(lm_biomass_un)


p1<-ggplot(noaa, aes(x = noaa$Calcifying, y = noaa$Herbivore, shape = noaa$INHABITATION, fill = noaa$INHABITATION))+
  geom_point()+
  stat_smooth(method = "lm")+
  theme_bw()+
  xlab("Calcifying benthic cover (%)")+
  ylab("Herbivore fish biomass (g/m2)")+
  scale_shape_manual(values = c(21, 24))+
  scale_fill_manual(values = c("Grey", "Black"))+
  guides(fill=FALSE, shape = FALSE)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  ggsave("Calcifying x Herbivore.pdf", width=8, height=8, units="cm")

p1

lm1b<- lm(noaa$Herbivore ~ noaa$HARD.CORAL)
summary(lm1b)

p1b<-ggplot(noaa, aes(x = noaa$HARD.CORAL, y = noaa$Herbivore, shape = noaa$INHABITATION, fill = noaa$INHABITATION))+
  geom_point()+
  stat_smooth(method = "lm")+
  theme_bw()+
  xlab("Coral cover (%)")+
  ylab("Herbivore fish biomass (g/m2)")+
  scale_shape_manual(values = c(21, 24))+
  scale_fill_manual(values = c("Grey", "Black"))+
  guides(fill=FALSE, shape = FALSE)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  ggsave("Coral x Herbivore.pdf", width=8, height=8, units="cm")

p1b

lm2<- lm(log10(noaa$Cells_per_ml) ~ noaa$Calcifying)
summary(lm2)

set.panel(2,2)
p2<-ggplot(noaa, aes(x = noaa$HARD.CORAL, y = log10(noaa$Cells_per_ml), shape = noaa$INHABITATION, fill = noaa$INHABITATION))+
  geom_point()+
  stat_smooth(method = "lm")+
  theme_bw()+
  xlab("Calcifying benthic cover (%)")+
  ylab("Log10 Cells per ml")+
  scale_shape_manual(values = c(21, 24))+
  scale_fill_manual(values = c("Grey", "Black"))+
  guides(fill=FALSE, shape = FALSE)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  ggsave("Coral x Cells per ml.pdf", width=8, height=8, units="cm")

p2

lm3<- lm(log10(noaa$Herbivore) ~ log10(noaa$Cells_per_ml))
summary(lm3)

p3<-ggplot(noaa, aes(x = log10(noaa$Cells_per_ml), y = log10(noaa$Herbivore), shape = noaa$INHABITATION, fill = noaa$INHABITATION))+
  geom_point()+
  stat_smooth(method = "lm")+
  theme_bw()+
  xlab("Log 10 Cells per ml")+
  ylab("log 10 Herbivore")+
  scale_shape_manual(values = c(21, 24))+
  scale_fill_manual(values = c("Grey", "Black"))+
  guides(fill=FALSE, shape = FALSE)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  ggsave("Herbivore x Cells per ml.pdf", width=8, height=8, units="cm")

p3

lm4<- lm(log10(noaa$Prey.Fish) ~ log10(noaa$Cells_per_ml))
summary(lm4)

p4<-ggplot(noaa, aes(x = log10(noaa$Cells_per_ml), y = log10(noaa$Prey.Fish), shape = noaa$INHABITATION, fill = noaa$INHABITATION))+
  geom_point()+
  stat_smooth(method = "lm")+
  theme_bw()+
  xlab("Log10 Cells per ml")+
  ylab("Prey fish (g/m2)")+
  scale_shape_manual(values = c(21, 24))+
  scale_fill_manual(values = c("Grey", "Black"))+
  #guides(fill=FALSE, shape = FALSE)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  ggsave("Prey fish x Cells per ml.pdf", width=8, height=8, units="cm")

p4

lm5<- lm(log10(noaa$Predator.fish) ~ log10(noaa$Prey.Fish))
summary(lm5)

p5<-ggplot(noaa, aes(x = log10(noaa$Prey.Fish), y = log10(noaa$Predator.fish), shape = noaa$INHABITATION, fill = noaa$INHABITATION))+
  geom_point()+
  stat_smooth(method = "lm")+
  theme_bw()+
  ylab("Predator fish (g/m2)")+
  xlab("Prey fish (g/m2)")+
  scale_shape_manual(values = c(21, 24))+
  scale_fill_manual(values = c("Grey", "Black"))+
  #guides(fill=FALSE, shape = FALSE)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  ggsave("Predator fish x Prey fish.pdf", width=8, height=8, units="cm")

p5

lm6<- lm(log10(noaa$VLP_per_ml_Corr) ~ log10(noaa$Cells_per_ml))
summary(lm6)

p6<-ggplot(noaa, aes(x = log10(noaa$Cells_per_ml), y = log10(noaa$VLP_per_ml_Corr), shape = noaa$INHABITATION, fill = noaa$INHABITATION))+
  geom_point()+
  stat_smooth(method = "lm")+
  theme_bw()+
  ylab("Log10 VLP per ml")+
  xlab("Log10 Cells per ml")+
  scale_shape_manual(values = c(21, 24))+
  scale_fill_manual(values = c("Grey", "Black"))+
  #guides(fill=FALSE, shape = FALSE)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  ggsave("VLP x Cells.pdf", width=8, height=8, units="cm")

p6

p7<-ggplot(noaa, aes(x = noaa$HARD.CORAL, y = log10(noaa$Herbivore), shape = noaa$INHABITATION, fill = noaa$INHABITATION))+
  geom_point()+
  stat_smooth(method = "lm")+
  theme_bw()+
  ylab("Log10 Herbivore")+
  xlab("Coral")+
  scale_shape_manual(values = c(21, 24))+
  scale_fill_manual(values = c("Grey", "Black"))+
  guides(fill=FALSE, shape = FALSE)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  ggsave("Herbivore x Coral.pdf", width=8, height=8, units="cm")

p7

p8<-ggplot(noaa, aes(x = noaa$HARD.CORAL, y = noaa$Shark, shape = noaa$INHABITATION, fill = noaa$INHABITATION))+
  geom_point()+
  stat_smooth(method = "lm")+
  theme_bw()+
  ylab("Shark")+
  xlab("Coral")+
  scale_shape_manual(values = c(21, 24))+
  scale_fill_manual(values = c("Grey", "Black"))+
  guides(fill=FALSE, shape = FALSE)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  ggsave("Shark x Coral.pdf", width=8, height=8, units="cm")

p8

p9<-ggplot(noaa, aes(x = noaa$HARD.CORAL, y = log10(noaa$VLP_per_ml_Corr), shape = noaa$INHABITATION, fill = noaa$INHABITATION))+
  geom_point()+
  stat_smooth(method = "lm")+
  theme_bw()+
  ylab("Log10 VLP")+
  xlab("Coral")+
  scale_shape_manual(values = c(21, 24))+
  scale_fill_manual(values = c("Grey", "Black"))+
  guides(fill=FALSE, shape = FALSE)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  ggsave("VLP x Coral.pdf", width=8, height=8, units="cm")

p9
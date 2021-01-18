library(ggplot2)

noaaov<-read.csv("Microbial_and_macrobial_variables.csv")
biomass<-read.csv("Biomass_data_with_inhabitation.csv")

#______
#MAP
library(mapdata)
library(ggplot2)
noaa.de<-read.csv("Microbial_and_macrobial_variables.csv", header = TRUE)
noaa.de$Cells_log10<-log10(noaa.de$Cells_per_ml)


mp1 <- fortify(map(fill=TRUE, plot=FALSE))
mp2 <- mp1
mp2$long <- mp2$long + 360
mp2$group <- mp2$group + max(mp2$group) + 1
mp <- rbind(mp1, mp2)

map<-ggplot(data = mp, aes(x = long, y = lat, group = group)) +
  geom_path()+   
  scale_x_continuous(limits = c(120, 280)) +
  scale_y_continuous(limits = c(-50, 50)) +
  theme_bw() +
  geom_point(data = noaa.de, aes(group=NULL,x = Degrees.East, y = LATITUDE, size=HARD.CORAL, color=Cells_log10))+
  scale_size_continuous(range = c(2,12))+
  scale_colour_gradient2(low = "deepskyblue2", mid = "grey", high = "darkred", midpoint = 6, space="Lab")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  ggsave("Map_new.eps", width = 26, height = 16, units = "cm" )
map

##################
#Separate maps for each variable

install.packages("viridis")
library("viridis")

jitter <- position_jitter(width = 3, height = 1)

#1 - Coral
map_coral<-ggplot(data = mp, aes(x = long, y = lat, group = group)) +
  geom_path()+   
  scale_x_continuous(limits = c(120, 280)) +
  scale_y_continuous(limits = c(-30, 30)) +
  theme_minimal() +
  geom_point(data = noaa.de, aes(group=NULL,x = Degrees.East, y = LATITUDE, fill=HARD.CORAL), shape = 21, size = 4, position=jitter)+
  scale_fill_viridis(option="magma", direction = -1) +
  #scale_fill_gradient2(low = "darkred", mid = "grey", high = "skyblue4", midpoint = 20, space="Lab")+
  theme(
    #panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    legend.title = element_blank(),
    legend.position=c(0.8, 0.5))+
  ggsave("Map_coral.eps", width = 26, height = 10, units = "cm" )
map_coral
#2 - Total fish biomass
map_fish<-ggplot(data = mp, aes(x = long, y = lat, group = group)) +
  geom_path()+   
  scale_x_continuous(limits = c(120, 280)) +
  scale_y_continuous(limits = c(-30, 30)) +
  theme_minimal() +
  geom_point(data = noaa.de, aes(group=NULL,x = Degrees.East, y = LATITUDE, fill=(Herbivore+Invertivore+Manta+Other+Planktivore+Piscivore+Shark)), shape = 21, size = 4, position=jitter)+
  scale_fill_viridis(option="magma", direction = -1) +
  #scale_fill_gradient2(low = "darkred", mid = "grey", high = "skyblue4", midpoint = 20, space="Lab")+
  theme(
    #panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    legend.title = element_blank(),
    legend.position=c(0.8, 0.5))+
  ggsave("Map_fish.eps", width = 26, height = 10, units = "cm" )
map_fish

#3 - Bacterial abundance
map_bac<-ggplot(data = mp, aes(x = long, y = lat, group = group)) +
  geom_path()+   
  scale_x_continuous(limits = c(120, 280)) +
  scale_y_continuous(limits = c(-30, 30)) +
  theme_minimal() +
  geom_point(data = noaa.de, aes(group=NULL,x = Degrees.East, y = LATITUDE, fill=(Cells_per_ml)), shape = 21, size = 4, position=jitter)+
  scale_fill_viridis(option="magma", direction = -1) +
  #scale_fill_gradient2(low = "darkred", mid = "grey", high = "skyblue4", midpoint = 20, space="Lab")+
  theme(
    #panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    legend.title = element_blank(),
    legend.position=c(0.8, 0.5))+
  ggsave("Map_bac.eps", width = 26, height = 10, units = "cm" )
map_bac


#3 - Bacterial abundance log10
map_bac_log<-ggplot(data = mp, aes(x = long, y = lat, group = group)) +
  geom_path()+   
  scale_x_continuous(limits = c(120, 280)) +
  scale_y_continuous(limits = c(-30, 30)) +
  theme_minimal() +
  geom_point(data = noaa.de, aes(group=NULL,x = Degrees.East, y = LATITUDE, fill=(log10(Cells_per_ml))), shape = 21, size = 4, position=jitter)+
  scale_fill_viridis(option="magma", direction = -1) +
  #scale_fill_gradient2(low = "darkred", mid = "grey", high = "skyblue4", midpoint = 20, space="Lab")+
  theme(
    #panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    legend.title = element_blank(),
    legend.position=c(0.8, 0.5))+
  ggsave("Map_bac_log.eps", width = 26, height = 10, units = "cm" )
map_bac_log


#______
#Overview of variable relationships

panel.cor <- function(x, y, digits=2, prefix="", cex.cor, ...) {
  usr <- par("usr")
  on.exit(par(usr))
  par(usr = c(0, 1, 0, 1))
  r <- abs(cor(x, y, use="complete.obs"))
  txt <- format(c(r, 0.123456789), digits=digits)[1]
  txt <- paste(prefix, txt, sep="")
  if(missing(cex.cor)) cex.cor <- 0.8/strwidth(txt)
  text(0.5, 0.5, txt, cex = cex.cor * (2 + r) / 2)
}

panel.hist <- function(x, ...) {
  usr <- par("usr")
  on.exit(par(usr))
  par(usr = c(usr[1:2], 0, 1.5) )
  h <- hist(x, plot = FALSE)
  breaks <- h$breaks
  nB <- length(breaks)
  y <- h$counts
  y <- y/max(y)
  rect(breaks[-nB], 0, breaks[-1], y, col="white", ...)
}


postscript(file="pairwise.plot.ps", width=14, height=12, horizontal = TRUE)
pairs(rfdataov, upper.panel = panel.cor,
      diag.panel = panel.hist,
      lower.panel = panel.smooth)
dev.off()

#Biomass 

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

#Other pairwise relationships NOT in the final manuscript

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

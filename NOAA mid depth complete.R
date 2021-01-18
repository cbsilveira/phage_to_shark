library(ggplot2)
library(fields)
library(Hmisc)


#_______________
#RandomForests
library(rfPermute)
library(randomForest)

noaaov<-read.csv("P2S_CRED_Data_OV.csv")
rfdataov<-noaaov[,c(15,16,18,19,22,23,27,28,31,32,33,35,36,37)]
rfdataov$VLP_per_ml_Corr<-log10(rfdataov$VLP_per_ml_Corr)
rfdataov$Cells_per_ml<-log10(rfdataov$Cells_per_ml)
rfdataov<-na.omit(rfdataov)

#Full dataset with oceanographic data - this has 110 sites 
set.seed(156)
summary(rfdataov)
pairs(rfdataov)
library(rpart)
rfdataov_tree<-rpart(rfdataov$HARD.CORAL~., data=rfdataov)
summary(rfdataov_tree)
plot(rfdataov_tree)
text(rfdataov_tree)


rftune<-tuneRF(rfdataov[,-4],rfdataov[,4])
rftune<-tuneRF(rfdataov[,-11],rfdataov[,11])

rf_allvars<-rfPermute(rfdataov$HARD.CORAL ~ ., data = rfdataov[,2:14], mtry = 8, proximity = TRUE, ntree = 1000, nrep = 1000)
rf_allvars
plot(rf_allvars)
plot(rp.importance(rf_allvars, scale = TRUE))
temp<- MDSplot(rf_allvars, noaa$REGION)

#remove the vars of low significance:
rfdataov_2nd<-rfPermute(rfdataov$HARD.CORAL ~ ., data = rfdataov[,c(2,3,4,5,6,7)], mtry = 4, proximity = TRUE, ntree = 1000, nrep = 1000)
rfdataov_2nd
plot(rfdataov_2nd)
plot(rp.importance(rfdataov_2nd, scale = TRUE))
temp_2nd<- MDSplot(rfdataov_2nd, noaa$REGION)


#further remove herbivore, at the bottom of the importance list
rfdataov_3<-rfPermute(rfdataov$HARD.CORAL ~ ., data = rfdataov[,c(2,3,4,5,6)], mtry = 4, proximity = TRUE, ntree = 1000, nrep = 1000)
rfdataov_3
plot(rfdataov_3)
plot(rp.importance(rfdataov_3, scale = TRUE))
temp_3<- MDSplot(rfdataov_3, noaa$REGION)

#Random forest with benthos only
rfdataov_4<-rfPermute(rfdataov$HARD.CORAL ~ ., data = rfdataov[,2:4], mtry = 2, proximity = TRUE, ntree = 1000, nrep = 1000)
rfdataov_4
plot(rfdataov_4)
plot(rp.importance(rfdataov_4, scale = TRUE))

#Random forest with no benthos
rfdataov_5<-rfPermute(rfdataov$HARD.CORAL ~ ., data = rfdataov[,5:14], mtry = 4, proximity = TRUE, ntree = 1000, nrep = 1000)
rfdataov_5
plot(rfdataov_5)
plot(rp.importance(rfdataov_5, scale = TRUE))
rfdataov_5$pval
rfdataov_5$importance


# No microbial data:
rf_nomic<-rfPermute(rfdataov$HARD.CORAL ~ ., data = rfdataov[,c(7:14)], mtry = 4, proximity = TRUE, ntree = 1000, nrep = 1000)
rf_nomic
plot(rf_nomic)
plot(rp.importance(rf_nomic, scale = TRUE))
temp_2nd<- MDSplot(rf_nomic, noaa$REGION)

#Separate Inhabited and Uninhabited

#1. Inhabited

inhabited_all<-noaaov[which(noaaov$INHABITATION == "Inhabited"),]
inhabited_rfdata_all<-inhabited_all[,c(15,16,18,19,22,23,27,28,31,32,33,35,36,37)]
inhabited_rfdata_all$VLP_per_ml_Corr<-log10(inhabited_rfdata_all$VLP_per_ml_Corr)
inhabited_rfdata_all$Cells_per_ml<-log10(inhabited_rfdata_all$Cells_per_ml)
irfd_all<-na.omit(inhabited_rfdata_all)

rf_inhabited<-rfPermute(irfd_all$HARD.CORAL ~ ., data = irfd_all[,c(5:14)], mtry = 4, proximity = TRUE, ntree = 1000, nrep = 1000)
rf_inhabited
plot(rp.importance(rf_inhabited, scale = TRUE))
rf_inhabited$pval
rf_inhabited$importance

#2. Uninhabited

uninhabited_all<-noaaov[which(noaaov$INHABITATION == "Uninhabited"),]
uninhabited_rfdata_all<-uninhabited_all[,c(15,16,18,19,22,23,27,28,31,32,33,35,36,37)]
uninhabited_rfdata_all$VLP_per_ml_Corr<-log10(uninhabited_rfdata_all$VLP_per_ml_Corr)
uninhabited_rfdata_all$Cells_per_ml<-log10(uninhabited_rfdata_all$Cells_per_ml)
urfd_all<-na.omit(uninhabited_rfdata_all)

rf_uninhabited<-rfPermute(urfd_all$HARD.CORAL ~ ., data = urfd_all[,c(5:14)], mtry = 4, proximity = TRUE, ntree = 1000, nrep = 1000)
rf_uninhabited
plot(rp.importance(rf_uninhabited, scale = TRUE))
rf_uninhabited$pval
rf_uninhabited$importance

###########################
#Conditional random forests

library(party)

crf<-cforest(rfdataov$HARD.CORAL ~ ., data = rfdataov[,2:14])
summary(crf)
crf_imp<-varimp(crf,  conditional = TRUE)
barplot(crf_imp)

crf_nobenthos<-cforest(rfdataov$HARD.CORAL ~ ., data = rfdataov[,5:14])
summary(crf_nobenthos)
crf_nobenthos_imp<-varimp(crf_nobenthos,  conditional = TRUE)
barplot(crf_nobenthos_imp)

############################

#Splines
#spline with benthic variables, and microbes (no herbivore)
tps_toprfvars<- Tps(rfdataov[,c(2,3,4,5,6)], rfdataov$HARD.CORAL)
summary(tps_toprfvars)
plot(tps_toprfvars)
surface(tps_toprfvars, type = "p")

#Plot Microbe and VLP:
out.p1<-predictSurface(tps_toprfvars, xy=c(4,5))
plot.surface(out.p1, type = "C",zlim = c(0,40))

#spline with VMR and fish (herbivores and predators)
rfdataov_nolog<-noaaov[,c(15,16,18,19,22,23,27,28,31,32,33,35,36,37)]
rfdataov_nolog<-na.omit(rfdataov_nolog)
rfdataov_nolog$VMR<-rfdataov_nolog$VLP_per_ml_Corr/rfdataov_nolog$Cells_per_ml
rfdataov_nolog$Predators<-rfdataov_nolog$Shark+rfdataov_nolog$Piscivore

tps_vmr_fish<- Tps(rfdataov_nolog[,c(7,15,16)], rfdataov_nolog$HARD.CORAL)
summary(tps_vmr_fish)
plot(tps_vmr_fish)
surface(tps_vmr_fish, type = "p")
out.tps_vmr_fish<-predictSurface(tps_vmr_fish, xy=c(2,3))
plot.surface(out.tps_vmr_fish, type = "C")


###########################
#Previous analyses with mtry = 3, Full dataset with oceanographic data - this has 100 sites because one datapoint is missing from one site
rfdov_2
rfdov_2$importance
rfdov_2$pval
set.panel(1,1)
plot(rfov_2)
summary(rfdov_2)
varImpPlot(rfov_2, type = 1)
plot(rp.importance(rfov_2, scale = TRUE))
set.panel(1,1)
temp<- MDSplot(rfov_2, noaa$REGION)

#No temp and DIC
rfdov_3<-rfPermute(rfdov$HARD.CORAL ~ ., data = rfdov[,c(2,3,4,5,6,7,8,9,10,11,13)], mtry = 3, proximity = TRUE, ntree = 1000, nrep = 1000)
rfdov_3
set.panel(1,1) 
varImpPlot(rfdov_3, type = 1)
plot(rp.importance(rfdov_3, scale = TRUE))

#No benthos, no temp and no DIC
rfdov_4<-rfPermute(rfdov$HARD.CORAL ~ ., data = rfdov[,c(5,6,7,8,9,10,11,13)], mtry = 3, proximity = TRUE, ntree = 1000, nrep = 1000)
rfdov_4
set.panel(1,1) 
varImpPlot(rfdov_4, type = 1)
plot(rp.importance(rfdov_4, scale = TRUE))

#No benthos
rfdov_5<-rfPermute(rfdov$HARD.CORAL ~ ., data = rfdov[,c(5,6,7,8,9,10,11,12,13,14)], mtry = 3, proximity = TRUE, ntree = 1000, nrep = 1000)
rfdov_5
set.panel(1,1) 
varImpPlot(rfdov_5, type = 1)
plot(rp.importance(rfdov_5, scale = TRUE))

#No oceanographic data, 111 sites
noaa<-read.csv("P2S_CRED Data.csv")

#data transformations (log10 cells and microbes)
rfdata<-noaa[,c(14,15,17,18,23,24,29,31,34,37,38)]
rfdata<-na.omit(rfdata)
rfdata$VLP_per_ml_Corr<-log10(rfdata$VLP_per_ml_Corr)
rfdata$Cells_per_ml<-log10(rfdata$Cells_per_ml)
rfd<-na.omit(rfdata)

set.seed(156)
rf1<- rfPermute(rfd$HARD.CORAL ~ ., data = rfd[,2:11], proximity = TRUE, ntree = 1000, nrep = 1000)
rf1
rf1$importance
rf1$pval
set.panel(1,1)
plot(rf1)
summary(rf1)
varImpPlot(rf1, type = 1)
proximityPlot(rf1, dim.x = 1, dim.y = 2, legend.loc = NULL, point.size = 2, circle.size = NULL,circle.border = NULL, hull.alpha = 0.3, plot = TRUE)
plot(rp.importance(rf1, scale = TRUE))
set.panel(1,1)
temp<- MDSplot(rf1, noaa$REGION)

# No microbial data:
rf1b<- rfPermute(rfd$HARD.CORAL ~ ., data = rfd[,c(2,3,4,7,8,9,10,11)], proximity = TRUE, ntree = 1000, nrep = 1000)
rf1b
plot(rf1b)
varImpPlot(rf1b, type = 1)
proximityPlot(rf1b, dim.x = 1, dim.y = 2, legend.loc = NULL, point.size = 2, circle.size = NULL,circle.border = NULL, hull.alpha = 0.3, plot = TRUE)
plot(rp.importance(rf1b, scale = TRUE))
set.panel(1,1)
temp<- MDSplot(rf1b, noaa$REGION)

#Because benthic variables are covariates, I removed them and tested again
rf4<- rfPermute(rfd$HARD.CORAL ~ ., data = rfd[,5:11], proximity = TRUE, ntree = 1000, nrep = 1000)
rf4
rf4$importance
rf4$pval
plot(rf4)
varImpPlot(rf4, type = 1)
proximityPlot(rf4, dim.x = 1, dim.y = 2, legend.loc = NULL, point.size = 2, circle.size = NULL,circle.border = NULL, hull.alpha = 0.3, plot = TRUE)
plot(rp.importance(rf4, scale = TRUE))
temp<- MDSplot(rf4, noaa$REGION)

# No microbial data:
rf4b<- rfPermute(rfd$HARD.CORAL ~ ., data = rfd[,c(7,8,9,10,11)], proximity = TRUE, ntree = 1000, nrep = 1000)
rf4b
set.panel(1,1)
plot(rf4b)
varImpPlot(rf4b, type = 1)
proximityPlot(rf4b, dim.x = 1, dim.y = 2, legend.loc = NULL, point.size = 2, circle.size = NULL,circle.border = NULL, hull.alpha = 0.3, plot = TRUE)
plot(rp.importance(rf4b, scale = TRUE))
temp<- MDSplot(rf4b, noaa$REGION)


#Separate Inhabited and Uninhabited

#1. Inhabited

inhabited<-noaa[which(noaa$INHABITATION == "Inhabited"),]
inhabited_rfddata<-inhabited[,c(14,15,17,18,23,24,29,31,34,37,38)]
inhabited_rfddata$VLP_per_ml_Corr<-log10(inhabited_rfddata$VLP_per_ml_Corr)
inhabited_rfddata$Cells_per_ml<-log10(inhabited_rfddata$Cells_per_ml)
irfd<-na.omit(inhabited_rfddata)

rf2<- rfPermute(irfd$HARD.CORAL ~ ., data = irfd[,5:11], proximity = TRUE, ntree = 1000, nrep = 1000)
rf2
rf2$importance
rf2$pval
plot(rf2)
summary(rf2)
varImpPlot(rf2, type = 1)
proximityPlot(rf2)
plot(rp.importance(rf2, scale = TRUE))

uninhabited<-noaa[which(noaa$INHABITATION == "Uninhabited"),]
uninhabited_rfddata<-uninhabited[,c(14,15,17,18,23,24,29,31,34,37,38)]
uninhabited_rfddata$VLP_per_ml_Corr<-log10(uninhabited_rfddata$VLP_per_ml_Corr)
uninhabited_rfddata$Cells_per_ml<-log10(uninhabited_rfddata$Cells_per_ml)
urfd<-na.omit(uninhabited_rfddata)

rf3<- rfPermute(urfd$HARD.CORAL ~ ., data = urfd[,5:11], proximity = TRUE, ntree = 1000, nrep = 1000)
rf3
rf3$importance
rf3$pval
plot(rf3)
summary(rf3)
varImpPlot(rf3, type = 1)
proximityPlot(rf3)
plot(rp.importance(rf3, scale = TRUE))

#__________
#Thin Plate Splines with cca, macroalage, turf, VLP, Cells, and Herbivores, no sharks 
spline_data<-rfd[,c(2,3,4,5,6,7)] #this includes benthos
tps1<- Tps(rfd[,c(2,3,4,5,6,7)], rfd$HARD.CORAL)
summary(tps1)
plot(tps1)
surface(tps1, type = "p")

rfd$HARD.CORAL
length(tps1$y)
rsquared_tps1<-(cor(rfd$HARD.CORAL,tps1$y))^2

rsquared_tps1


#Plot Microbe and Herbivores:
out.p1<-predictSurface(tps1, xy=c(5,6))
plot.surface(out.p1, type = "C", zlim = c(0,100))

#Plot Microbe and VLP:
set.panel(1,1)
out.p2<-predictSurface(tps1, xy=c(5,4))
plot.surface(out.p2, type = "C", zlim = c(0,100))

#Plot Microbes and herbivores 
out.p2s<-predictSurface(tps1, xy=c(5,6))
plot.surface(out.p2s, type = "C", zlim = c(0,100))

#Plot TA and Microbes:
out.p3<-predictSurface(tps1, xy=c(3,5))
plot.surface(out.p3, type = "C", zlim = c(0,100))

#Plot MA and Microbes:
out.p4<-predictSurface(tps1, xy=c(2,5))
plot.surface(out.p4, type = "C", zlim = c(0,100))

#Plot CCA and Herbivore:
out.p5<-predictSurface(tps1, xy=c(1,6))
plot.surface(out.p5, type = "C", zlim = c(0,100))

#Plot VMP and Hervivore
out.p6<-predictSurface(tps1, xy=c(4,6))
plot.surface(out.p6, type = "C", zlim = c(0,100))

set.panel(2,2)
#Microbes x VLP
plot.surface(out.p2, type = "C",zlim = c(0,100))
#Microbes x Herbivore
plot.surface(out.p1, type = "C",zlim = c(0,100))
#Microbes x Turf algae
plot.surface(out.p3, type = "C",zlim = c(0,100))
#Microbes x Macroalgae
plot.surface(out.p4, type = "C",zlim = c(0,100))


set.panel(3,1)
#Microbes x VLP
plot.surface(out.p2, type = "C",zlim = c(0,100))
#VLP x Herbivore
plot.surface(out.p6, type = "C", zlim = c(0,100))
#Cells x Herbivore
plot.surface(out.p1, type = "C",zlim = c(0,100))

##
#Thin Plate Splines with macroalage, turf, VLP, Cells, Herbivores, and predators (piscivores and sharks) - the removal of CCA is arbitrary, because it is the benthic variable with least importance
spline_data_nocca<-rfd_large[,c(3,4,5,6,7,12)] #this includes benthos
tps1b_nocca<- Tps(rfd_large[,c(3,4,5,6,7,12)], rfd$HARD.CORAL)
summary(tps1b_nocca)
plot(tps1b_nocca)
surface(tps1b_nocca, type = "p")

#Plot Herbivores and predators:
out.tps1b_nocca<-predictSurface(tps1b_nocca, xy=c(5,6))
plot.surface(out.tps1b_nocca, type = "C", zlim = c(0,100))
#Plot Cells and VLP:
out.tps1b_nocca_2<-predictSurface(tps1b_nocca, xy=c(4,3))
plot.surface(out.tps1b_nocca_2, type = "C", zlim = c(0,100))
#plot VLP and herbivore
out.tps1b_nocca_3<-predictSurface(tps1b_nocca, xy=c(5,3))
plot.surface(out.tps1b_nocca_3, type = "C", zlim = c(0,100))


########

#Thin Plate Splines with macroalgae, turf algae, herbivores, predators (piscivores+sharks), and VMR(VLP/Cells)

rfd_large<-rfd
rfd_large$predators<-rfd$Shark+rfd$Piscivore
rfd_large$VMR<-rfdata$VLP_per_ml_Corr/rfdata$Cells_per_ml
tps2<- Tps(rfd_large[,c(3,4,7,12,13)], rfd_large$HARD.CORAL)
summary(tps2)
tps2
plot(tps2)
surface(tps2, type = "p")

#Herbivore and predator
out.tps2_1<-predictSurface(tps2, xy=c(3,4))
plot.surface(out.tps2_1, type = "C")

#VMR and predator fish
out.tps2_2<-predictSurface(tps2, xy=c(4,5))
plot.surface(out.tps2_2, type = "C")

###
#Thin Plate Splines with VLP, Cells, Herbivores, predators, alkalinity, and planktivore (no dic and temp)
bio_alk_data<-rfdov[,c(1,5,6,7,8,9,10,11,13)]
bio_alk_data$predators<-bio_alk_data$Piscivore+bio_alk_data$Shark
tps5<- Tps(bio_alk_data[,c(2,3,4,6,9,10)], bio_alk_data$HARD.CORAL)
summary(tps5)
plot(tps5)
surface(tps5, type = "p")

#Plot Microbe and VLPs:
surface1<-predictSurface(tps5, xy=c(1,2))
plot.surface(surface1, type = "C", zlim = c(0,100))

#Plot VLP and predators:
surface2<-predictSurface(tps5, xy=c(1,6))
plot.surface(surface2, type = "C", zlim = c(0,100))

#Plot Microbes and predators:
surface3<-predictSurface(tps5, xy=c(2,6))
plot.surface(surface3, type = "C", zlim = c(0,100))

#Plot Microbes and alkalinity:
surface4<-predictSurface(tps5, xy=c(2,5))
plot.surface(surface4, type = "C", zlim = c(0,100))

###
#Thin Plate Splines with VLP, Cells, Herbivores, predators, alkalinity, DIC and planktivore (no temp)
sig_data_nobenthos<-rfdov[,c(1,5,6,7,8,9,10,11,12,13)]
sig_data_nobenthos$predators<-sig_data_nobenthos$Piscivore+sig_data_nobenthos$Shark
tps6<- Tps(sig_data_nobenthos[,c(2,3,4,6,9,11)], sig_data_nobenthos$HARD.CORAL) #I removed alkalinity, the least important in this list, because tps only takes 6 vars
summary(tps6)
plot(tps6)
surface(tps6, type = "p")

#Plot Microbe and VLPs:
surface6.1<-predictSurface(tps6, xy=c(1,2))
plot.surface(surface6.1, type = "C", zlim = c(0,100))

#Plot VLP and predators:
surface6.2<-predictSurface(tps6, xy=c(1,6))
plot.surface(surface6.2, type = "C", zlim = c(0,100))

#Plot Microbes and predators:
surface6.3<-predictSurface(tps6, xy=c(2,6))
plot.surface(surface6.3, type = "C", zlim = c(0,100))

#Plot Microbes and alkalinity:
surface6.4<-predictSurface(tps6, xy=c(2,5))
plot.surface(surface6.4, type = "C", zlim = c(0,100))


#########

#Robust and Quantile smoothing using a thin-plate spline
qtps1<- QTps(rfd[,c(2,3,4,5,6,7)], rfd$HARD.CORAL)
set.panel(2,2)
plot(qtps1)
set.panel(2,2)
surface(qtps1, type = "p")
qtps1$fitted.values
summary(qtps1)
plot(qtps1)
qtps1$residuals
qtps1$fitted.values

#smoothing by quantiles:
qtps1_q2<- QTps(rfd[,c(2,3,4,5,6,7)], rfd$HARD.CORAL,alpha =0.2)
summary(qtps1_q2)

qtps1_q4<- QTps(rfd[,c(2,3,4,5,6,7)], rfd$HARD.CORAL,alpha =0.4)
summary(qtps1_q4)

qtps1_q6<- QTps(rfd[,c(2,3,4,5,6,7)], rfd$HARD.CORAL,alpha =0.6)
summary(qtps1_q6)

qtps1_q6<- QTps(rfd[,c(2,3,4,5,6,7)], rfd$HARD.CORAL,alpha =0.8)
summary(qtps1_q8)

qtps1_q100<- QTps(rfd[,c(2,3,4,5,6,7)], rfd$HARD.CORAL,alpha =1)
summary(qtps1_q10)

#Plot effective DFs for each quantile
quantile_splines<- data.frame("Quantiles" = c(20,40,60,80,100),"Eff.df" = c(qtps1_q2$eff.df, qtps1_q4$eff.df, qtps1_q6$eff.df, qtps1_q8$eff.df,qtps1_q10$eff.df))
fit.smspl_quantiles = smooth.spline(quantile_splines$Eff.df,quantile_splines$Quantiles, cv=TRUE)

quant_dfs<-ggplot(quantile_splines, aes( x = Quantiles, y = Eff.df))+
  geom_point(shape = 21, size = 4, stroke = 1)+
  stat_smooth(method = "glm", formula = y ~ poly(x, 4), se = FALSE)+
  #lines(fit.smspl_quantiles, col="red")+
  #geom_line(linetype = "dashed")+
  theme_classic()+
  ylim(75,115)+
  ylab("Effective Degrees of Freedom")+
  xlab("Calcifying cover quantiles")+
  guides(fill=FALSE)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  ggsave("Eff_df_quantiles.eps", width = 9, height = 9, units = "cm" )
quant_dfs  

coral_cover_boxplot<-ggplot(noaa.de, aes(y = METHOD, x = HARD.CORAL, fill = HARD.CORAL))+
  geom_boxplot(outlier.shape = NA, width = 0.7)+
  geom_point(position=position_jitter(width=0, height=0.4),  shape = 21)+
  scale_fill_viridis(option="magma", direction = -1)+
  theme_classic()+
  theme(
    legend.position = "none",
    #axis.title.y=element_blank(),
    axis.title.x=element_blank(),
    axis.text.y = element_blank(),
    axis.text.x = element_text(size=rel(0.9)))+
  ggsave("Coral_cover.eps", width = 8.5, height = 3, units = "cm" )
coral_cover_boxplot


#Run quantile regression splines with each variable separatly and obtain the DFs to see which ones are contributing to bumpiness:
#Microbes:
library(fields)
qsreg1<-qsreg(rfd$Cells_per_ml, rfd$HARD.CORAL)
qsreg1$alpha
summary(qsreg1)
set.panel(2,2)
plot(qsreg1)
#VLP:
qsreg2<-qsreg(rfd$VLP_per_ml_Corr, rfd$HARD.CORAL)
summary(qsreg2)
plot(qsreg2)
#Herbivore:
qsreg3<-qsreg(rfd$Herbivore, rfd$HARD.CORAL)
summary(qsreg3)
plot(qsreg3)
#Macroalgae:
qsreg4<-qsreg(rfd$MA, rfd$HARD.CORAL)
summary(qsreg4)
plot(qsreg4)
#Predator:
qsreg5<-qsreg(rfd$Shark, rfd$HARD.CORAL)
summary(qsreg5)
plot(qsreg5)
#PredatorxPrey:
qsreg6<-qsreg(rfd$Shark, rfd$Herbivore)
summary(qsreg6)
plot(qsreg6)
#TA:
qsreg8<-qsreg(rfd$TA, rfd$HARD.CORAL)
summary(qsreg8)
plot(qsreg8)
#CCA:
qsreg9<-qsreg(rfd$CCA, rfd$HARD.CORAL)
summary(qsreg9)
plot(qsreg9)

qsreg7<-qsreg(rfd$VLP_per_ml_Corr, rfd$Cells_per_ml)
summary(qsreg7)
plot(qsreg7)

##Smoothing spline code from Toni

fit.smspl1 = smooth.spline(rfd$Cells_per_ml, rfd$HARD.CORAL,cv=TRUE)
fit.smspl2 = smooth.spline(rfd$VLP_per_ml_Corr, rfd$HARD.CORAL,cv=TRUE)
fit.smspl3 = smooth.spline(rfd$Herbivore, rfd$HARD.CORAL, cv = TRUE)
fit.smspl4 = smooth.spline(rfd$Shark, rfd$HARD.CORAL,cv=TRUE)

set.panel(2,2)
plot(rfd$Cells_per_ml, rfd$HARD.CORAL)+lines(fit.smspl1, col="red")
plot(rfd$VLP_per_ml_Corr, rfd$HARD.CORAL)+lines(fit.smspl2, col="red")
plot(rfd$Herbivore, rfd$HARD.CORAL)+lines(fit.smspl3, col="red")
plot(rfd$Shark, rfd$HARD.CORAL)+lines(fit.smspl4, col="red")

#BY FISHING

fit.smspl1u = smooth.spline(urfd$Cells_per_ml, urfd$HARD.CORAL,cv=TRUE)
fit.smspl1u
plot(urfd$HARD.CORAL~urfd$Cells_per_ml)+lines(fit.smspl1u, col="red")

fit.smspl1i = smooth.spline(irfd$Cells_per_ml,irfd$HARD.CORAL,cv=TRUE)
fit.smspl1i
plot(irfd$HARD.CORAL~irfd$Cells_per_ml)+lines(fit.smspl1i, col="red")


fit.smspl5 = smooth.spline(urfd$Shark, urfd$HARD.CORAL, cv = TRUE)
plot(urfd$HARD.CORAL~urfd$Shark)+lines(fit.smspl5, col="red")

fit.smspl6 = smooth.spline(irfd$Shark, irfd$HARD.CORAL,cv=TRUE)
plot(irfd$HARD.CORAL~irfd$Shark)+lines(fit.smspl6, col="red")

#________
#Multiple Linear regression

mlm1<- lm(HARD.CORAL ~ Shark + Cells_per_ml + VLP_per_ml_Corr + Herbivore + TA + MA + CCA, data = noaa)
summary(mlm1)
plot(mlm1)
anova(mlm1)

#Calculate relative importance of variables
library(relaimpo)
calc.relimp(mlm1,type=c("lmg","last","first","pratt"),rela=TRUE)

# Bootstrap Measures of Relative Importance (1000 samples) 
boot <- boot.relimp(mlm1, b = 1000, type = c("lmg", 
                                             "last", "first", "pratt"), rank = TRUE, 
                    diff = TRUE, rela = TRUE)
booteval.relimp(boot) # print result
plot(booteval.relimp(boot,sort=TRUE)) # plot result

#Step the lm to find best model

library(MASS)
step <- stepAIC(mlm1, direction="both")
step$anova

mlm2<- lm(HARD.CORAL ~ Cells_per_ml + VLP_per_ml_Corr + Herbivore + TA + 
            MA + CCA, data = noaa)
summary(mlm2)

set.panel(1,1)
plot.surface(out.p1, type = "p")

library(plot3D)
set.panel(1,1)
scatter3D(log10(noaa$VLP_per_ml_Corr),log10(noaa$Cells_per_ml), noaa$HARD.CORAL) 

#________
#GAM - GLM
library(mgcv)
gam_y <- gam( HARD.CORAL ~ s(Shark, bs = 'cr') + s(Cells_per_ml, bs = 'cr') + s(VLP_per_ml_Corr,bs = 'cr') + s(Herbivore,bs = 'cr' ) + s(TA,bs = 'cs') + s(MA,bs = 'cr') + s(CCA,bs = 'cr'), data = rfd, method = "REML")
summary(gam_y)
plot(gam_y)
par(mfrow = c(2,2))
gam.check(gam_y)

#Calculate relative importance of variables
library(relaimpo)
calc.relimp(gam_y,type=c("lmg","last","first","pratt"),rela=TRUE)

# Bootstrap Measures of Relative Importance (1000 samples) 
boot_gam <- boot.relimp(gam_y, b = 1000, type = c("lmg", 
                                             "last", "first", "pratt"), rank = TRUE, 
                    diff = TRUE, rela = TRUE)
booteval.relimp(boot_gam) # print result
plot(booteval.relimp(boot_gam,sort=TRUE)) # plot result

#Step the lm to find best model
library(gam)
library(MASS)
step_gam <- step.Gam(gam_y, direction="both")
step_gam$anova

#gam with vars from RF no benthos
gam_yb <- gam( HARD.CORAL ~ s(Shark, bs = 'cr') + s(Cells_per_ml, bs = 'cr') + s(VLP_per_ml_Corr,bs = 'cr') + s(Herbivore,bs = 'cr'), data = rfd, method = "REML")
summary(gam_yb)
plot(gam_yb)
par(mfrow = c(2,2))
gam.check(gam_yb)
# Bootstrap Measures of Relative Importance (1000 samples) 
boot_gamyb <- boot.relimp(gam_yb, b = 1000, type = c("lmg", 
                                                  "last", "first", "pratt"), rank = TRUE, 
                        diff = TRUE, rela = TRUE)
booteval.relimp(boot_gamyb) # print result
plot(booteval.relimp(boot_gamyb,sort=TRUE)) # plot result

#Explanations for methods

#lmg is the R2 contribution averaged over orderings among regressors, cf. e.g. Lindeman, Merenda and Gold 1980, p.119ff or Chevan and Sutherland (1991).
#last is each variables contribution when included last, also sometimes called usefulness.
#first is each variables contribution when included first, which is just the squared covariance between y and the variable.
#pratt is the product of the standardized coefficient and the correlation

#______
#MAP
library(mapdata)
library(ggplot2)
noaa.de<-read.csv("P2S_CRED Data_degrees_east.csv", header = TRUE)
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
#I decided to make separate maps for each variable

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

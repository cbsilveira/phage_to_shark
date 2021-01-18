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

#Random forest with no benthos
rfdataov_5<-rfPermute(rfdataov$HARD.CORAL ~ ., data = rfdataov[,5:14], mtry = 4, proximity = TRUE, ntree = 1000, nrep = 1000)
rfdataov_5
plot(rfdataov_5)
plot(rp.importance(rfdataov_5, scale = TRUE))
rfdataov_5$pval
rfdataov_5$importance

#Random forest with benthos only
rfdataov_4<-rfPermute(rfdataov$HARD.CORAL ~ ., data = rfdataov[,2:4], mtry = 2, proximity = TRUE, ntree = 1000, nrep = 1000)
rfdataov_4
plot(rfdataov_4)
plot(rp.importance(rfdataov_4, scale = TRUE))

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



#########

#Robust and Quantile smoothing using a thin-plate spline
qtps1<- QTps(rfdataov[,c(2,3,4,5,6,7)], rfdataov$HARD.CORAL)
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
qtps1_q2<- QTps(rfdataov[,c(2,3,4,5,6,7)], rfdataov$HARD.CORAL,alpha =0.2)
summary(qtps1_q2)

qtps1_q4<- QTps(rfdataov[,c(2,3,4,5,6,7)], rfdataov$HARD.CORAL,alpha =0.4)
summary(qtps1_q4)

qtps1_q6<- QTps(rfdataov[,c(2,3,4,5,6,7)], rfdataov$HARD.CORAL,alpha =0.6)
summary(qtps1_q6)

qtps1_q6<- QTps(rfdataov[,c(2,3,4,5,6,7)], rfdataov$HARD.CORAL,alpha =0.8)
summary(qtps1_q8)

qtps1_q100<- QTps(rfdataov[,c(2,3,4,5,6,7)], rfdataov$HARD.CORAL,alpha =1)
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
qsreg1<-qsreg(rfdataov$Cells_per_ml, rfdataov$HARD.CORAL)
qsreg1$alpha
summary(qsreg1)
set.panel(2,2)
plot(qsreg1)
#VLP:
qsreg2<-qsreg(rfdataov$VLP_per_ml_Corr, rfdataov$HARD.CORAL)
summary(qsreg2)
plot(qsreg2)
#Herbivore:
qsreg3<-qsreg(rfdataov$Herbivore, rfdataov$HARD.CORAL)
summary(qsreg3)
plot(qsreg3)
#Macroalgae:
qsreg4<-qsreg(rfdataov$MA, rfdataov$HARD.CORAL)
summary(qsreg4)
plot(qsreg4)
#Predator:
qsreg5<-qsreg(rfdataov$Shark, rfdataov$HARD.CORAL)
summary(qsreg5)
plot(qsreg5)
#PredatorxPrey:
qsreg6<-qsreg(rfdataov$Shark, rfdataov$Herbivore)
summary(qsreg6)
plot(qsreg6)
#TA:
qsreg8<-qsreg(rfdataov$TA, rfdataov$HARD.CORAL)
summary(qsreg8)
plot(qsreg8)
#CCA:
qsreg9<-qsreg(rfdataov$CCA, rfdataov$HARD.CORAL)
summary(qsreg9)
plot(qsreg9)

qsreg7<-qsreg(rfdataov$VLP_per_ml_Corr, rfdataov$Cells_per_ml)
summary(qsreg7)
plot(qsreg7)




#________
#Multiple Linear regression

mlm1<- lm(HARD.CORAL ~ Shark + Cells_per_ml + VLP_per_ml_Corr + Herbivore + TA + MA + CCA, data = noaaov)
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
scatter3D(log10(noaaov$VLP_per_ml_Corr),log10(noaaov$Cells_per_ml), noaaov$HARD.CORAL) 

#________
#GAM - GLM
library(mgcv)
gam_y <- gam( HARD.CORAL ~ s(Shark, bs = 'cr') + s(Cells_per_ml, bs = 'cr') + s(VLP_per_ml_Corr,bs = 'cr') + s(Herbivore,bs = 'cr' ) + s(TA,bs = 'cs') + s(MA,bs = 'cr') + s(CCA,bs = 'cr'), data = rfdataov, method = "REML")
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
gam_yb <- gam( HARD.CORAL ~ s(Shark, bs = 'cr') + s(Cells_per_ml, bs = 'cr') + s(VLP_per_ml_Corr,bs = 'cr') + s(Herbivore,bs = 'cr'), data = rfdataov, method = "REML")
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

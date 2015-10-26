rm(list=ls(all=TRUE))
library(dplyr); library(survival)
detach("package:plyr", unload=TRUE) ## This line was added to prevent getting errors when you run n()
set.seed(1) ## to make sure parametric boostrap gives same values every time
setwd('C:/Users/Joseph/Documents/SempaetalAJE-00426-2015')
patients <- read.csv(file='patients.csv', header = T)

obs<-patients[,c('id','end.yrs','start.yrs','end.int','start.int')] ## creating a file from perdat.edmb with only 3 columns
#####################################################
#I've trimmed out all the covariates and only simulate viral set point, and
#a viral load that varies around the set point and affects the actual risk.
#########################################################
mod.p<-as.matrix(NA, nr=1000, nc=7)
mod.pvalue<-as.data.frame(mod.p)
colnames(mod.pvalue)[1]<-'mod1.cf'
mod.pvalue$mod1.pv<-NA
mod.pvalue$mod2.cf<-NA
mod.pvalue$mod2.pv<-NA
mod.pvalue$mod3.cf<-NA
mod.pvalue$mod3.pv<-NA
mod.pvalue$mod.events<-NA
numSims <- 1000

for(ii in 1:numSims){
threshVL <- 400
#set.seed(2017)
#browser()

obs <- mutate(group_by(obs, id)
              , setPoint = rlnorm(1, meanlog=7, sdlog=1)
              , VL = rlnorm(n(), meanlog=log(setPoint), sdlog=1)
              , VL = pmax(VL, threshVL)
              , logVL = log(VL)
              , cumVL = cumsum(logVL)
)

summary(obs)

##################################################################
#Calculate hazards following the Cox assumptions, as far as I can tell.
##################################################################

base <- 0.1
bet_VL <- 0.2
bet_cum <- 0

obs <- mutate(obs
              , timeHazard = base*(end.yrs-start.yrs)
              , lp = bet_VL*logVL + bet_cum*cumVL
              , hazard = timeHazard*exp(lp)
              , OI = rbinom(n(), size=1, prob=1-exp(-hazard))
)
summary(obs)
# rdsave(obs)
###################################################################
#Survival model
#Make a survival model, and try to recapture the real dynamics 
#(logVL contributes linearly to the linear predictor, and cumVL should 
#have no residual effect). The model just gets it wrong.
#########################################################################
coxph(Surv(start.yrs, end.yrs, OI) ~ logVL, data=obs)
mod1<-coxph(Surv(start.yrs, end.yrs, OI) ~ cumVL + logVL, data=obs)
mod.cont<- coef(summary(mod1))
sum.mod1<-summary(mod1)
mod.events<-as.numeric(sum.mod1$nevent)
mat.cont<-as.matrix(mod.cont[,1])
mod1.cf<-as.numeric(mat.cont[1,1])
mat.cont<-as.matrix(mod.cont[,5])
mod1.pv<-as.numeric(mat.cont[1,1])
#########################################################################
#glm, on the other hand, nails it.
#########################################################################
obs <- within(obs, {
  dur <- end.yrs - start.yrs
  ldur <- log(dur)
})

summary(glm(
  OI ~ cumVL + logVL,
  offset = ldur,
  family=binomial(link=cloglog),
  data=obs
))
mod2<-glm(formula = OI ~ cumVL + logVL, family = binomial(link = cloglog), 
    data = obs, offset = ldur)
mod.glm<- coef(summary(mod2))
mat.glm<-as.matrix(mod.glm[,1])
mod2.cf<-as.numeric(mat.glm[2,1])
mat.glm<-as.matrix(mod.glm[,4])
mod2.pv<-as.numeric(mat.glm[2,1])
##########################################################################
#Note: Rounding to eighths instead gets rid of most of the missingness, but also lets cumVL 
#start creeping backing in.
############################################################################
#########################################################################
coxph(Surv(start.int, end.int, OI) ~ logVL, data=obs)
mod3<-coxph(Surv(start.int, end.int, OI) ~ cumVL + logVL, data=obs)
mod.int<- coef(summary(mod3))
mat.int<-as.matrix(mod.int[,1])
mod3.cf<-as.numeric(mat.int[1,1])
mat.int<-as.matrix(mod.int[,5])
mod3.pv<-as.numeric(mat.int[1,1])
mod.pvalue[ii,]<-rbind(matrix(c(mod1.cf,mod1.pv,mod2.cf,mod2.pv,mod3.cf,mod3.pv,mod.events),nrow=1,ncol=7))
}

#########################################################################
#Ploting histograms of p-values from the different models
#########################################################################

for(ftype in c('pdf', 'eps')) {
# tiff(file.path('Figures','distribution of p-Values for cVL generated from 3 models.tiff'),
#      w = 480, h = 1020, units="px", pointsize=18)
# labels <- c('CoxPH continuous interval', 'Poisson reg continuous interval', 'CoxPH rounded interval')
    breaks <- seq(0,1, by = .025)
    mains <- c('A)', 'B)', 'C)')
    if(ftype=='pdf') { ## if pdf
        pdf(paste0('Figures/distribution of p-Values for cVL generated from 3 models.pdf'), w = 8.5, h = 8.5, pointsize = 10)
        par(mfrow=c(3,1))
    }
    
  for(ii in 1:3) {
    if(ftype=='eps') { ##creating figure files in eps format accoridng to AJE requirements
      postscript(file.path('Figures',paste0('AJE-00426-2015 Sempa Web Fig 5',ii,'.eps')), onefile = F, width = 8.5, height = 8.5, horizontal = F, pointsize = 12)
      par(mfrow=c(1,1), mar = c(4,4.5,1.5,1.5), mgp = c(3,1,0), oma = c(1,0,0,0), 'ps'=12)
    }
#     hist(mod.pvalue[[paste0('mod',ii,'.pv')]], xlim = c(0,1), main = labels[ii], col = 'black', breaks = breaks,
#            xlab = 'P value for cumVL effect (when there is none)')
    hist(mod.pvalue[[paste0('mod',ii,'.pv')]], xlim = c(0,1), col = 'black', breaks = breaks, main = '', xlab = '')
      abline(v=0.05,untf=F, col="blue", lty=4, lwd=4) ## line showing statistical significance
      if(ii==1) mtext(expression(P~'-'~value~of~cVL[2]~from~Cox~PH~model~with~unrounded~intervals), side = 1, line = 2)
      if(ii==2) mtext(expression(P~'-'~value~of~cVL[2]~from~Poisson~model~with~unrounded~intervals), side = 1, line = 2)
      if(ii==3) mtext(expression(P~'-'~value~of~cVL[2]~from~Cox~PH~model~with~rounded~intervals), side = 1, line = 2)
      if(ftype=='pdf') mtext(mains[ii], 3, line = 0, adj=0)
        if(ftype=='eps') dev.off() ## if eps end figure inside loop
  }
  
  if(ftype!='eps')  dev.off() ## if pdf end figure
}

summary((mod.pvalue$mod.events))

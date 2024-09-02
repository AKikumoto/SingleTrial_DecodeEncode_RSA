# Conjunction_Simulate
# simulating regression parts (data generated from RSA models)
# TO DO
# - Extend this to do dimensioality simulation??? (change model function?)
# --> test relationship between conjunction and dimensionality!!
# - think of axies and make this indexing flexible
# =============================================================
#Removing evrything from workspace
graphics.off()
rm(list = ls(all = TRUE))

#Setting up directory
# Dir_M<-dirname(file.choose())# pick whichever file of RSA.txt model
Dir_R <- path.expand("~/Dropbox/w_ONGOINGRFILES/w_OTHERS")
Dir_M <- path.expand("~/Dropbox/w_SCRIPTS/P.HAL2017/HAL2017_BI/DataAnalysis/R/w_Simulation")

# Load libraries
library(data.table)
library(tidyverse)
library(broom)
library(foreach)
library(doMC)
library(doParallel)
library(scales)
library(forcats)
library(stringr)
library(psych)
library(lazyeval)
library(caret)

# Source Files
setwd(Dir_R)
source('basic_lib.R')
setwd(Dir_M)
source('Conj_Sim_Envs.R')

# --------------------------------
# Simulation settings
# --------------------------------
# Possible things to modulate (in loop)
sess <- "expand"
niter <- 500;# number of iteractions
noise_u <- seq(0.5,3,0.1);# amount of noise
# nobs <- seq(1200,1600,200);# number of observations
nobs <- 2000; # average (Exp1 = 2300, Exp2 = 2071.875)

# Summarize
axisALL <- expand.grid(noise_u,nobs)
names(axisALL) <- c("noise_u","nobs")

# --------------------------------
# Load RSA Models and base units
# --------------------------------
g(mL, uL)  %=%  loadModel_Units(file.path(Dir_M,"CONJ_v2"))
# g(mL, uL) %=% preprop_decode(ds_bl,ds_eeg,prepSet)

# --------------------------------
# Loop for simulation 
# --------------------------------
registerDoMC(detectCores()-2);
# ni<-1;ii<-1

# Nested parallelized loop
results<-foreach(ni=1:niter,.combine='cbind')%:%
  foreach(ii=1:dim(axisALL)[1],.combine='c') %dopar% {
    # Get settings (make this part flexible...)
    noise_u <- axisALL[["noise_u"]][ii]
    nobs <- axisALL[["nobs"]][ii]
    
    # --------------------------------
    # Generate traininng/test data
    # --------------------------------
    genL <-list(
      ncontext = nrow(uL[[1]]), # number of classes
      ntrial = nobs, # number of observations
      uL = uL, # list of base units defined for each feature in RSA
      nnoises = nrow(mL[[1]]), # number of non-informative units
      noiseF_g = function(n) {rnorm(n,0,1.5)}, # noise in noisy units
      noiseF_u = function(n) {rnorm(n,0,noise_u)},  # noise in informative units
      onoff = c("RULE","STIM","RESP") # name of features to include (need to match names in mL)
    )

    ds <- genData(genL)

    # --------------------------------
    # Run classification & Summarize results
    # --------------------------------
    control<-trainControl(method="repeatedcv",
                          number=4,# number of folds
                          repeats=8,# number of repetitions
                          p=0.75,# proportion of data in the fold?
                          selectionFunction = "oneSE",
                          classProbs=T,allowParallel=T,savePredictions=T)


    f <- as.formula('LABEL ~ .') # should be always same
    m <- train(f, data=ds, method="pda", metric="Accuracy", trControl=control);
    g(acc, pred) %=% sort_decode(ds,m,varL) # rehaped pred data frame

    # --------------------------------
    # RSA
    # --------------------------------
    modelL <-list(
      pred = pred,#prediction from classification model
      mL = mL, #list of all possible RSA models
      logit = T,#whether if we convert probability into logit
      rsa = names(mL) #list of RSA moedls to activate: Ex) c("CONJ","RULE","STIM","RESP")
    )

    r <- modelRSAFit(modelL)

    # --------------------------------
    # Summarize
    # --------------------------------
    r_s <- describeBy(r$tvalue,r$vars)
    r_s <- bind_rows(lapply(r_s,as.data.table),.id="term") %>% data.table()
    r_s$acc <- acc$Accuracy # decoding accuracy
    r_s[,nsim:=ni] # loop index of number of iterations
    r_s[,sidx:=ii] # loop index of simulation settings
    r_s[,biasB:=(mean > sd)*1] # bias bigger than 1SD
    return(list(r_s))
}

# Get all results
rr <- bind_rows(results)
rAGG <- rr[,lapply(.SD,mean),.SDcols=c("mean","sd","acc","biasB"),by=c("sidx","term")]
rAGG <- cbind(rAGG,axisALL[rAGG$sidx,]);print(rAGG)
saveRDS(rAGG, paste0(sess,"_simulatedDimScore.rds"))

# Plot (heatmap)
# Heat map: time-frequencuy results
quartz(width=6.5,height=4.5);theme_set(theme_bw(base_size = 20))#32/28
ggplot() +
  geom_tile(data=rAGG[term=="CONJ_RSA"],aes(x=nobs,y=noise_u,fill=biasB),linejoin = "round")+
  #geom_contour(data=r_accG,aes(time,freqL,z=DV),colour='white',alpha=1,bins=1)
  #Annotation -------------------
  ylab("Noise")+xlab("# of trials")+labs(fill="DV")+
  #Axis + scales-------------------
  scale_x_continuous(expand=c(0,0))+
  #scale_y_discrete(expand=c(0,0),breaks=c(70,90,120,150))+
  #scale_y_discrete(expand=c(0,0),breaks=c(4,12,20,50))+
  scale_fill_gradientn(limits=c(-0.02,1),colours=viridis::inferno(30),oob=squish,n.breaks=4)+
  theme(panel.border=element_blank(),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        axis.line=element_line(size=0.75),
        legend.key = element_blank(),
        legend.key.size=unit(1,"lines"),
        legend.position = "right",
        legend.title=element_blank(),
        legend.text =element_text(family="Helvetica",size=16),
        axis.ticks.length = unit(0.5, "lines"),
        axis.title  = element_text(family="Helvetica",vjust=1.8),
        plot.title = element_text(hjust = 0, size = 10),
        strip.background=element_blank())

# # Plot (histogram)
# quartz(width=5,height=6)
# ggplot(r, aes(x=tvalue,group=vars)) +
#   geom_vline(xintercept = 0)+
#   facet_wrap(~vars) +
#   #coord_cartesian(xlim=c(-10,10),ylim=c(1,500))+
#   xlab("t-values")+
#   geom_histogram(binwidth=0.2,color="black", fill="white")


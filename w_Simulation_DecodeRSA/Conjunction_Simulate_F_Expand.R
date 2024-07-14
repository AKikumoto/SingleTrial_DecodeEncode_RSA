# Conjunction_Simulate expanded space
# simulating regression parts (data generated from RSA models)
# =============================================================
#Removing evrything from workspace
graphics.off()
rm(list = ls(all = TRUE))

#Setting up directory
# Dir_RSAMODEL<-dirname(file.choose())# pick whichever file of RSA.txt model
Dir_HERE<-dirname(rstudioapi::getSourceEditorContext()$path)
Dir_RSAMODEL<-file.path(Dir_HERE,"CONJ_v2")

# Load libraries
library(data.table)
library(tidyverse)
library(broom)
library(foreach)
library(doMC)
library(scales)
library(forcats)
library(stringr)
library(psych)
library(caret)

# Classification setting
#"pda" = penalized linear discriminant analysis (what we actually used)
#"svmRadial" = support vector machine
#"nb" = naive bayes
#"rf" = random forest
#"svmLinear3" = L2 Regularized Support Vector Machine (dual) with Linear Kernel
method<-'pda';
formula<-as.formula(paste('LABEL ~ .'))
metric<-"Accuracy";
control<-trainControl(method="repeatedcv",
                      number=2,repeats=8,p=0.75,
                      selectionFunction = "oneSE",
                      classProbs=TRUE,allowParallel=TRUE,savePredictions=TRUE)

# Local functions
is.nanM <- function(x)do.call(cbind, lapply(x, is.nan))

#======================================================================================================
#Load RSA Models (adjusted for different decoding variables)
#======================================================================================================
setwd(Dir_RSAMODEL)

# Models 
mlist <- list.files(pattern = ".txt")
mL<-lapply(as.list(mlist),function(f){read.table(f,header = F)})
for (m in 1:length(mlist)){assign(gsub(".txt","",mlist)[m],mL[[m]])}

# Hypothetical units
mL_S<-lapply(mL,function(x){x<-x/x;x[is.nanM(x)]<-0;x});# put all units in a same range of activation
names(mL_S)<-c("CONJ_H","CONJ_L","ODDEVEN","RESP","RULE","STIM")
for (v in names(mL_S)){names(mL_S[[v]])<-print(paste0(v,"_",names(mL_S[[v]])))}

#======================================================================================================
# Generate data using hypothetical units (12 units for each feature = 48 units)
#======================================================================================================
registerDoMC(6);i<-1
ntrial<-2000
ncontext<-dim(CONJ_H_M)[1];
ds<-data.table(LABEL=rep(1:ncontext,1,ntrial))

r<-foreach(i=1:ntrial) %dopar%{
  # 12 hypothetical units (V1~V12) for each feature
  cl <- ds$LABEL[i] #cl <- sample(ds$LABEL,1);# permutation test
  rule_s <- mL_S[["RULE"]][cl,] + rnorm(ncontext,0,0.5); 
  stim_s <- mL_S[["STIM"]][cl,] + rnorm(ncontext,0,0.5);
  resp_s <- mL_S[["RESP"]][cl,] + rnorm(ncontext,0,0.5);
  oddeven_s <- mL_S[["ODDEVEN"]][cl,] + rnorm(ncontext,0,0.5);
  conj_low_s <- mL_S[["CONJ_M"]][cl,] + rnorm(ncontext,0,0.5);
  conj_high_s <- mL_S[["CONJ_H_M"]][cl,] + rnorm(ncontext,0,0.5);
  noise <- data.table(t(rnorm(ncontext,0,1.5)));
  #noise <- data.table(t(runif(ncontext,min=0,max=1)));
  names(noise) <- paste0("NOISE_",names(noise)) 
  obs <- cbind(rule_s,stim_s,resp_s,oddeven_s,noise); # signal (Adjust this part!)
}

ds <- cbind(ds,bind_rows(r))
ds[, LABEL:=LETTERS[ds$LABEL]]

#======================================================================================================
# #Run classification
#======================================================================================================
# When classifications are perfect, the rest will break so be careful!
m <- train(formula, data=ds, method=method, metric=metric, trControl=control);

# #Summarize results
# pred <- summary_decode(m,varL,list(PRED=T,CM=F,IMP=F,ROC=F))$r_prob
pred<-m$pred %>% data.table();
pred[,acc:=(pred==obs)*1]
pred[,c("pred","obs","lambda","Resample"):=NULL]
pred<-pred[,lapply(.SD, mean),by="rowIndex"][order(rowIndex)]

# #Reshape data and prepar
mesVs <- colnames(pred)[colnames(pred) %in% LETTERS]
idVs <- colnames(pred)[!colnames(pred) %in% c(mesVs,"acc")]
castF1 <- paste0("CLASS~",paste0(idVs[!idVs %in% c("SUBID")],collapse="+"))
castF2 <- paste0(paste0(idVs,collapse="+"),"~vars")
pred <- data.table::melt(pred,id.vars = idVs,measure.vars=mesVs,variable.name="CLASS",value.name="PP")
pred <- merge(rowid_to_column(ds,'rowIndex')[,c("rowIndex","LABEL")],pred,by="rowIndex")

# # Logit transform
pred[,obs:=as.numeric(as.factor(pred$LABEL))]
pred[PP==0,PP:=1e-10]
pred[PP==1,PP:=1-(1e-10)]
pred[,PP:=psych::logit(PP)]

#======================================================================================================
# RSA
#======================================================================================================
# Simulate n times
cc<-1
clist <-unique(pred$obs)

results<-
  foreach(cc=clist) %dopar% {
    d<-pred[obs==cc,];
    predF<-data.table::dcast(d,castF1, value.var="PP")
    
    # Random model
    #RAND_M <- CONJ_M[sample.int(nrow(CONJ_M)),]
    
    # This part could be more abstract with "eval" 
    # TSCONJ (Rule x Resp/Stim) with CJ RT control
    predF$RULE_RSA<-RULE_M[,cc];#RULE MODEL
    predF$STIM_RSA<-STIM_M[,cc];#STIM MODEL
    predF$RESP_RSA<-RESP_M[,cc];#RESP MODEL
    predF$ODDEVEN_RSA<-ODDEVEN_M[,cc];#ODDEVEN MODEL
    predF$CONJ_RSA<-CONJ_M[,cc];#CONJ_M MODEL
    predF$CONJ_H_RSA<-CONJ_H_M[,cc];#CONJ_M_H_M MODEL
    #predF$CTR_RSA<-RAND_M[,cc];#RANDOM MODEL
    dvIdx<-grepl("CLASS|*._RSA$",names(predF))
    
    # Regression all trials at once (but individually)!
    Y<-as.matrix(predF[,which(!dvIdx),with=FALSE]);#DVs(this indexing takes time...)
    X<-as.matrix(predF[,which(dvIdx)[-1],with=FALSE]);#IVs(first var is CLASS!)
    # r<-coeff(.lm.fit(cbind(1,X),Y));#fastest, but only gives unstandardized coefficients
    r<-ls.print(lsfit(X,Y),print.it=F);#includes intercept!
    
    # Summarize results!
    names(r$coef.table)<-colnames(predF)[!dvIdx]
    estList<-rownames(r$coef.table[[1]]);#list of vars for estimates
    r<-rbindlist(lapply(r$coef.table,as.data.frame),idcol=TRUE)
    setnames(r,old=c(".id","t-value","Pr(>|t|)","Std.Err"),new=c("rowIndex","tvalue","pvalue","SE"))
    r[,vars:=rep(estList,dim(r)[1]/length(estList))]
    r<-r[vars!="Intercept",c("rowIndex","vars","tvalue"),with=F]
    
    return(r)
  }

# Summarize
rr<-data.table(do.call(rbind,results))

# Summary stats by group
describeBy(rr$tvalue,rr$vars)

# Plot
quartz(width=5,height=6)
ggplot(rr, aes(x=tvalue,group=vars)) + 
  geom_vline(xintercept = 0)+
  facet_wrap(~vars) + 
  #coord_cartesian(xlim=c(-10,10),ylim=c(1,500))+
  xlab("t-values")+
  geom_histogram(binwidth=0.2,color="black", fill="white")


# Conjunction_Simulation_Environment
# simulating regression parts (data generated from RSA models)
# =============================================================

# Local functions
is.nanM <- function(x)do.call(cbind, lapply(x, is.nan))


# =============================================================
# Loading RSA models
# INPUT: dir = directory that contaitns the list of RSA models
# =============================================================
loadModel_Units <-function(dir){
  # Go to specified locations and read up txt files
  setwd(dir)
  mlist <- list.files(pattern = ".txt")
  mL<-lapply(as.list(mlist),function(f){read.table(f,header = F)})
  names(mL)<-gsub("_M.txt","",mlist)
  
  # List of RSA models
  for (m in 1:length(mlist)){assign(gsub(".txt","",mlist)[m],mL[[m]],envir=.GlobalEnv)}
  
  # List of base units (0 or 1 activation)
  uL<-lapply(mL,function(x){x<-x/x;x[is.nanM(x)]<-0;x});
  names(uL)<-gsub("_M.txt","",mlist)
  #for (v in names(mL_S)){names(mL_S[[v]])<-paste0(v,"_",names(mL_S[[v]]))}
  return(list(mL,uL))
}

# =============================================================
# Generating data with basic (RSA derived units)
# =============================================================
genData <- function(dsL){
  # Generate template
  ds <- data.table(LABEL=rep(1:dsL$ncontext,1,dsL$ntrial))
  
  # Informative units + noise
  for (f in dsL$onoff){
    units <- dsL$uL[[f]][ds$LABEL,]
    noise <- dsL$noiseF_u(dsL$ncontext*dsL$ntrial)
    dim(noise) <- c(dsL$ntrial,dsL$ncontext)
    units <- units + noise
    names(units) <- gsub("V",paste0(f,"_U"),names(units),)
    ds <- cbind(ds,units)
  }
  
  # General noise units
  noise <- dsL$noiseF_g(dsL$ncontext*dsL$ntrial);
  dim(noise) <- c(dsL$ntrial,dsL$ncontext)
  noise <- data.table(noise)
  names(noise) <- gsub("V",paste0("NOISE_U"),names(noise),)
  ds <- cbind(ds,noise)
  
  # Relabel class label 
  ds[, LABEL:=LETTERS[ds$LABEL]]
  return(ds)
}

# =============================================================
# Summarize model output
# =============================================================
sort_decode <- function(ds,m,varL) {
  #Overall accuracy
  d<-m$trainingData
  best_pda<-function(m){m$results[m$results$lambda==as.character(m$bestTune),]};
  r_acc<-best_pda(m) %>% data.table();
  
  #Posterior probability (be careful to use data.table... rowIndex?)
  r_prob<-m$pred %>% data.table();
  r_prob[,acc:=(r_prob$pred==r_prob$obs)*1]# keep simple accuracy
  r_prob[,c("pred","obs","lambda","Resample"):=NULL]
  r_prob<-r_prob[,lapply(.SD, mean),by="rowIndex"][order(rowIndex)]
  
  # Reshape 
  mesVs <- colnames(r_prob)[colnames(r_prob) %in% LETTERS]
  idVs <- colnames(r_prob)[!colnames(r_prob) %in% c(mesVs,"acc")]
  r_prob <- melt(r_prob,id.vars=idVs,measure.vars=mesVs,variable.name="CLASS",value.name="PP")
  r_prob <- merge(rowid_to_column(ds,'rowIndex')[,c("rowIndex","LABEL")],r_prob,by="rowIndex")
  r_prob[,obs:=as.numeric(as.factor(r_prob$LABEL))]
  
  #Final outputs
  return(list(r_acc=r_acc,r_prob=r_prob))
}

# =============================================================
# Fitting model 
# =============================================================

modelRSAFit <- function(modelL){
  # Put outside of list for convenction
  pred <- modelL$pred
  #pred[,obs:=as.numeric(as.factor(pred$LABEL))]

  # Logit transform?
  if (modelL$logit){
    pred[PP==0,PP:=1e-10]
    pred[PP==1,PP:=1-(1e-10)]
    pred[,PP:=psych::logit(PP)]
  }

  # Do regression!
  r <- foreach(cc=unique(pred$obs)) %do% { #CAUTION! dopar does now work with lazy_eval method
      # Add predictors
      predF<-data.table::dcast(pred[obs==cc,],"CLASS~rowIndex",value.var="PP")
      if (any(grepl("RANDOM",modelL$rsa))){RAND_M <- CONJ_M[sample.int(nrow(CONJ_M)),]}
      for (m in modelL$rsa){predF[,(paste0(m,"_RSA")):=modelL$mL[[m]][,cc]]}
      #for (m in modelL$rsa){predF[,(paste0(m,"_RSA")):=lazy_eval(paste0(m,"_M[,cc]"),.GlobalEnv)]}
      dvIdx<-grepl("CLASS|*._RSA$",names(predF))

      # Regression all trials at once (but individually)!
      # r<-coeff(.lm.fit(cbind(1,X),Y));#fastest, but only gives unstandardized coefficients
      Y<-as.matrix(predF[,which(!dvIdx),with=FALSE]);#DVs(this indexing takes time...)
      X<-as.matrix(predF[,which(dvIdx)[-1],with=FALSE]);#IVs(first var is CLASS!)
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
  r <- data.table(do.call(rbind,r))
  return(r)
}

# =============================================================
# Others
# =============================================================

## %=% parsing list output to multiple variables 
# %-%, extendToMatch(), g()
# EX) g(a, b, c)  %=%  list("hello", 123, list("apples, oranges"))

# Generic form
'%=%' = function(l, r, ...) UseMethod('%=%')

# Binary Operator
'%=%.lbunch' = function(l, r, ...) {
  Envir = as.environment(-1)
  
  if (length(r) > length(l))
    warning("RHS has more args than LHS. Only first", length(l), "used.")
  
  if (length(l) > length(r))  {
    warning("LHS has more args than RHS. RHS will be repeated.")
    r <- extendToMatch(r, l)
  }
  
  for (II in 1:length(l)) {
    do.call('<-', list(l[[II]], r[[II]]), envir=Envir)
  }
}

# Used if LHS is larger than RHS
extendToMatch <- function(source, destin) {
  s <- length(source)
  d <- length(destin)
  
  # Assume that destin is a length when it is a single number and source is not
  if(d==1 && s>1 && !is.null(as.numeric(destin)))
    d <- destin
  
  dif <- d - s
  if (dif > 0) {
    source <- rep(source, ceiling(d/s))[1:d]
  }
  return (source)
}

# Grouping the left hand side
g = function(...) {
  List = as.list(substitute(list(...)))[-1L]
  class(List) = 'lbunch'
  return(List)
}

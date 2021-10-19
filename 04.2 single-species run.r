
### --------------------------------------------------------------------------------- ###
### ----------------------------------- Setup --------------------------------------- ###
### --------------------------------------------------------------------------------- ###

results_path <- "results/single-species/"

library(FLCore) # (version: ??2.5.20160107??)
library(FLash)  # (version: ??2.5.2??)
library(FLAssess)  # (version: ??2.5.20130716??)
library(FLXSA) # (version: ??2.5.20140808??)
library(ggplotFL) # (version: ??2.5.20160119??)
library(snow) # for multicore use (version: ??0.4.1??)
# library(forecast)


### --------------------------------------------------------------------------------- ###
### --------------------------------- load data ------------------------------------- ###
### --------------------------------------------------------------------------------- ###

load("data/cod-sim(1988).Rdata")        # assessed stock
load("data/survey idx.Rdata") # survey index
load("R/data/SRR.Rdata")      # SRR
source("R/functions/functions_incrM.R")


### --------------------------------------------------------------------------------- ###
### ------------------------------ General parameters ------------------------------- ###
### --------------------------------------------------------------------------------- ###

# assessCV <- 0.38 # Error in assessing stock numbers
# idxCV <- 0.25
idxCV <- 0

nyear <- 50
niter <- 1
year_seq <- seq(range(stkReal)["maxyear"]+1, length.out = nyear)
wts.nyears=5
fbar.nyears=5
ncores <- 1
# ncores <- as.numeric(Sys.getenv("NUMBER_OF_PROCESSORS"))
nn <- split(1:niter, 1:niter)

### --------------------------------------------------------------------------------- ###
### ------------------------------ Initial Assessment ------------------------------- ###
### --------------------------------------------------------------------------------- ###


# Conduct initial XSA
xsaCtrl <- FLXSA.control(tol = 1e-09, maxit = 80, min.nse = 0.3, fse = 1,
                         rage = 0, qage = 6, shk.n = TRUE, shk.f = TRUE,
                         shk.yrs = 5, shk.ages= 6, window = 100, tsrange = 99,
                         tspower = 1)


xsaReal <- FLXSA(stkReal, idxReal, xsaCtrl)
stkReal <- stkReal + xsaReal

# Record index catchability ====
q.hat <- xsaReal@q.hat 

# determine SRR (pay attention to the time window 88+ and 98+)
stkSr.name <- list(c("1988"),c("1998"))


### --------------------------------------------------------------------------------- ###
### ---------------------------- Create Objects for MSE ----------------------------- ###
### --------------------------------------------------------------------------------- ###

# Expand objects and duplicate observed objects ====
stkReal <- stf(stkReal, nyear=nyear, wts.nyears=wts.nyears, fbar.nyears=fbar.nyears) # used stf to carry forward biol indicies
idxReal <- window(idxReal, end=tail(year_seq,1))

# expand to niter
stkReal <- propagate(stkReal, niter)
for(i in seq(idxReal)){
  idxReal[[i]] <- propagate(idxReal[[i]], niter)
}
dims(idxReal[[1]])

# duplicate observed
stkObs <- stkReal
idxObs <- idxReal

# Create TAC object ====
TAC <- catch(stkReal) # TAC


### --------------------------------------------------------------------------------- ###
### ----------------------------------- M2 dynamics --------------------------------- ###
### --------------------------------------------------------------------------------- ###

m1 <- 0.2

###  1.  Seal(No.15) as predator
load("R/data/m2_seal.Rdata")
attach(m2_seal)
seal_1 <- mean(annual_PM2[Year %in% (2009:2013) & Prey.age == 1])  # 0.1310056 on age 1
seal_2 <- mean(annual_PM2[Year %in% (2009:2013) & Prey.age == 2])  # 0.1310056 on age 2
detach(m2_seal)

###  2. porpoise(No.16) as predator
load("R/data/m2_porpoise.Rdata")
attach(m2_porpoise)
porpoise_1 <- mean(annual_PM2[Year %in% (2009:2013) & Prey.age == 1])  # 0.4360056 on age 1
porpoise_2 <- mean(annual_PM2[Year %in% (2009:2013) & Prey.age == 2])  # 0.4360056 on age 2
detach(m2_porpoise)

###  3. cannibalism
fit_coef <- read.csv("R/data/fit_cod_coef.csv")

# in this case the max age is 6
fit_coef <- fit_coef[c(1:5,10:12),]
fit_coef$Group <- factor(fit_coef$Group)

# only age 1 and 2 as prey, m2 for age three maintains as that of 2014

##  4  Whiting(No.19) as predator
load("R/data/m2_whiting.Rdata")
attach(m2_whiting)
whiting_1 <- mean(annual_PM2[Year %in% (2009:2013) & Prey.age == 1])  # 0.006586642 on age 1
detach(m2_whiting)


### --------------------------------------------------------------------------------- ###
### --------------------------------- Stochasticity  -------------------------------- ###
### --------------------------------------------------------------------------------- ###

# Create errors for assessment and recruitment ====
set.seed(1)
recrErr <- vector(mode="list", length(stkSr))
for(i in seq(recrErr)){
  recrErr[[i]] <- FLQuant(exp(sample(c(residuals(stkSr[[i]])), nyear*niter, replace=TRUE)),
                          dimnames = list(year=year_seq, iter=1:niter))
}

idxErr <- idxReal
for(i in seq(idxReal)){
  idxErr[[i]]@index <- FLQuant(rep(1, prod(dim(idxReal[[i]]@index))), dimnames = dimnames(idxReal[[i]]@index))
  idxErrDim <- dim(idxReal[[i]]@index[,ac(year_seq)])
  idxErr[[i]]@index[,ac(year_seq)] <- rlnorm(prod(idxErrDim), meanlog = 0, sdlog = idxCV)
}
sd(log(c(idxErr[[1]]@index[3,ac(year_seq),,,,]))); idxCV # a check


### --------------------------------------------------------------------------------- ###
### ------------------------------------- MSE --------------------------------------- ###
### --------------------------------------------------------------------------------- ###


### -----------------   Building scenarios
loop_multi <- expand.grid(F=seq(0.2,0.8,0.1), srr=stkSr.name, cannibalism=c("canni","no canni"), update_m=c("update","no update"))
loop_multi <- loop_multi[1:42,]
loop_multi$scenarios <- paste(loop_multi$srr,loop_multi$F,loop_multi$cannibalism,loop_multi$update_m)


### ------------------- Loop to MSE

rn=41;x=1;ii=1

# for (rn in seq(nrow(loop_multi))){
for (rn in c(13,41)){
    
  srr = which(stkSr.name %in% loop_multi[rn,"srr"])
  F   = loop_multi[rn,"F"]
  
  # make duplicates of objects
  stkReal_dup  <- stkReal
  stkObs_dup   <- stkObs
  idxObs_dup   <- idxObs
  idxReal_dup  <- idxReal
  idxErr_dup   <- idxErr
  TAC_dup      <- TAC
  
  # with cannibalism
  if (loop_multi[rn,"cannibalism"] == "canni" ) {
    m(stkReal)[1,ac(year_seq),,,,] <- m1 + seal_1 + porpoise_1 + whiting_1
    m(stkReal)[2,ac(year_seq),,,,] <- m1 + seal_2 + porpoise_2 
  }
  
  # MSE function to run for each iter ====
  MSE_fun <- function(x){
    # load packages and functions for each cluster
    library(FLCore)
    library(FLash) 
    library(FLAssess) 
    library(FLXSA)
    
    # create iter versions of objects
    stkObs_nn <- iter(stkObs, x)
    stkReal_nn <- iter(stkReal, x)
    idxObs_nn <- iterIndicies(idxObs, x)
    idxReal_nn <- iterIndicies(idxReal, x)
    TAC_nn <- iter(TAC, x)
    idxErr_nn <- iterIndicies(idxErr, x)
    # stkErr_nn <- iter(stkErr, x)
    
    ii=ii+1
    for(ii in seq(year_seq)){
      # Assess stock from last years catch and indices
      if (ii > 1) {
        stkObs_nn <- assessStock(
          stock = stkObs_nn, index = idxObs_nn,
          control = xsaCtrl,
          maxyear = year_seq[ii]-1,
          method="XSA"
        )
      }
      
      # advance operational model by one year using TAC (save temporarily)
      tmp.n <-  advanceStock(
        stock = stkReal_nn, tac = TAC_nn, sr = stkSr[[srr]], 
        maxyear = year_seq[ii]-1,
        sr.residuals = iter(trim(recrErr[[srr]], year=year_seq[ii]), x),
        sr.residuals.mult = FALSE)
      
      # update m2 by cannibalism
      if (loop_multi[rn,"cannibalism"] == "canni") {
        # i=1;a=1
        for (i in 1:length(fit_coef$Group)) {           # extract every model
          for (a in 1:6) {                              # fill in m(stkReal) by age (row)
            
            if (i %in% 1:5 & a<=1) {                    # control the accumulative pm2 for prey age 1, predator age is i+1
              m(stkReal_nn)[a,ac(year_seq[ii]),,,,] = m(stkReal_nn)[a,ac(year_seq[ii]),,,,] + 
                exp(fit_coef[i,2] + 
                      fit_coef[i,3] * log(stock.n(tmp.n)[i+1,ac(year_seq[ii]),,,,]) + 
                      fit_coef[i,4] * log(stock.n(tmp.n)[a,ac(year_seq[ii]),,,,]))
            }
            
            if (i %in% 6:8 & a<=2) {                    # control the accumulative pm2 for prey age 2, predator age is i-2
              m(stkReal_nn)[a,ac(year_seq[ii]),,,,] = m(stkReal_nn)[a,ac(year_seq[ii]),,,,] + 
                exp(fit_coef[i,2] + 
                      fit_coef[i,3] * log(stock.n(tmp.n)[i-2,ac(year_seq[ii]),,,,]) + 
                      fit_coef[i,4] * log(stock.n(tmp.n)[a,ac(year_seq[ii]),,,,]))
            }
          }
        }
      }
      
      # set the M for next year temporarily for TAC (except for the last year)
      # if (ii < nyear)  {
      #   m(stkReal_nn)[,ac((year_seq[ii]+1)),,,,] <- m(stkReal_nn)[,ac((year_seq[ii])),,,,]
      # }
      
      # record the m into stkObs if we will have perfect knowledge of predation
      if (loop_multi[rn,"update_m"] == "update") {
        if (ii < nyear){m(stkObs_nn)[,ac(year_seq[ii]:(year_seq[ii]+1)),,,,] <-  m(stkReal_nn)[,ac(year_seq[ii]),,,,]}
      }
      
      # Calculate TAC if year ii is not the final simulation year
      if(ii < nyear){
        TAC_nn <- calcTac(stock=stkObs_nn, tac=TAC_nn, hcrFun=hcrConstF, yr_ass = year_seq[ii]-1, harvest=F,
                          stockSub_sr = stkSr[[srr]])
      }
      
      # advance operational model by one year using TAC (for stkReal)
      stkReal_nn <- advanceStock(
        stock = stkReal_nn, tac = TAC_nn, sr = stkSr[[srr]],
        maxyear = year_seq[ii]-1,
        sr.residuals = iter(trim(recrErr[[srr]], year=year_seq[ii]), x),
        sr.residuals.mult = FALSE
      )
      
      # record observed catches
      stkObs_nn <- observeCatch(stockReal=stkReal_nn, stockObs=stkObs_nn, obsErr = NULL)
      
      # record real and observed indices
      idxReal_nn <- observeIndex(stock=stkReal_nn, index=idxReal_nn, q.hat=q.hat, obsErr=NULL, minyear = year_seq[1])
      idxObs_nn <- observeIndex(stock=stkReal_nn, index=idxObs_nn, q.hat=q.hat, obsErr=idxErr_nn, minyear = year_seq[1])
      
      print(paste("iter =", x, "; year =", year_seq[ii]))
    }
    
    stkObs_nn <- assessStock(
      stock = stkObs_nn, index = idxObs_nn,
      control = xsaCtrl,
      maxyear = year_seq[ii],
      method="XSA"
    )
    
    RES <- list(
      stkObs_nn = stkObs_nn,
      stkReal_nn = stkReal_nn,
      idxObs_nn = idxObs_nn,
      idxReal_nn = idxReal_nn,
      TAC_nn = TAC_nn
    )
    
    sink(file="output.txt", append = TRUE)
    print(paste("run", rn, "iter", x, "completed @", Sys.time()))
    sink()
    
    return(RES)
  }  # end MSE
  
  #############        MSE running and result exportatin        ################
  cl <- makeCluster(ncores, type="SOCK", outfile="tmp.txt")
  clusterExport(cl, list=ls(), envir = .GlobalEnv)
  
  t1 <- Sys.time()
  res <- parLapply(cl, nn, MSE_fun )
  stopCluster(cl)
  Sys.time()- t1
  
  # Rebuild objects with multicore results (by iter) ====
  for(n in nn){
    iter(stkReal, n) <- res[[n]]$stkReal_nn
    iter(stkObs, n) <- res[[n]]$stkObs_nn
    iter(TAC, n) <- res[[n]]$TAC_nn
    # for(lev in seq(res[[n]]$idxReal_nn)){
    #   iter(idxReal[[lev]], n) <- res[[n]]$idxReal_nn[[lev]]
    #   iter(idxObs[[lev]], n) <- res[[n]]$idxObs_nn[[lev]]
    # }
  }
  
  
  
  ### Output / Analysis -------------------------------------------------------
  
  # Save / Load Results ====
  
  save( stkObs,  file=paste(results_path, "stkObs", loop_multi[rn,"scenarios"],".Rdata", sep="_") ) 
  save( stkReal, file=paste(results_path, "stkReal",loop_multi[rn,"scenarios"],".Rdata", sep="_") )
  # save( idxObs,  file=paste(results_path, "idxObs", loop_multi[rn,"scenarios"],".Rdata", sep="_") )
  # save( idxReal, file=paste(results_path, "idxReal",loop_multi[rn,"scenarios"],".Rdata", sep="_") )
  save( TAC,     file=paste(results_path, "TAC",    loop_multi[rn,"scenarios"],".Rdata", sep="_") )
  
  # restore objects
  stkReal  <- stkReal_dup
  stkObs   <- stkObs_dup
  idxObs   <- idxObs_dup
  idxReal  <- idxReal_dup
  idxErr   <- idxErr_dup
  TAC      <- TAC_dup
  
}



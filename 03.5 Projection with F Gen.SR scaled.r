
library(FLCore)
library(FLash)
library(FLBRP)
library(FLXSA)
library(FLAssess)
library(snow)

load(file = "R/data/gom.xsa.R")
load(file = "R/data/survey.R")
load(file = "R/data/SRR.xsa.R")      # SR
load(file = "R/data/q.hat.R")        # q.hat

source("R/functions/functions for GoM cod.R")
gom <- gom.xsa
ssb(gom)[,ac(2014),,,,]
stock.n(gom)[,ac(2014),,,,] <- stock.n(gom)[,ac(2014),,,,]*1449/1883.8

### --------------------------------- General Setting ---------------------------------- ###
### ------------------------------------------------------------------------------------ ###

###--------------------------- Fundemental Parameters ------------------------------ ###

niter <- 200
nyear <- 6
ncores <- 4 #as.numeric(Sys.getenv("NUMBER_OF_PROCESSORS"))
nn <- split(1:niter, 1:niter)
idxCV <- 0.12 

year_seq <- seq(range(gom)["maxyear"]+1, length.out = (nyear))

###--------------------------- Input data Processing ------------------------------ ###

# Expand objects and duplicate observed objects used stf to carry forward biol indicies
mtf.xsa <- stf(gom, nyear=nyear, wts.nyears=3, fbar.nyears=3) 
indices <- window(indices, end=tail(year_seq,1))

# expand to niter
mtf.xsa <- propagate(mtf.xsa, niter)
for(i in seq(indices)){
  indices[[i]] <- propagate(indices[[i]], niter)
}

# control object of XSA
xsa.control <- FLXSA.control(tol = 1e-09, maxit = 100, min.nse = 0.3, fse = 1,
                             rage = 0, qage = 8, shk.n = TRUE, shk.f = FALSE,
                             shk.yrs = 5, shk.ages= 5, window = 100, tsrange = 99,
                             tspower = 1)

# make the corresponding observational data available for assessment
mtf.xsa.obs <- mtf.xsa
indices.obs <- indices
TAC <- catch(mtf.xsa)

###--------------------------- Stochasticity ------------------------------ ###

set.seed(1)

# error for recruitment
SR.Err <- FLQuant(exp(sample(c(residuals(SR.xsa)[,ac(2002:2014),,,,]), nyear*niter, replace=TRUE)),
                  dimnames = list(year=year_seq, iter=1:niter))

# error for survey indices
indices.Err <- indices
for(i in seq(indices.Err)){
  indices.Err[[i]]@index <- FLQuant(rep(1, prod(dim(indices[[i]]@index))), dimnames = dimnames(indices[[i]]@index))
  indices.Err.Dim <- dim(indices[[i]]@index[,ac(year_seq)])
  indices.Err[[i]]@index[,ac(year_seq)] <- rlnorm(prod(indices.Err.Dim), meanlog = 0, sdlog = idxCV)
}
sd(log(c(indices.Err[[1]]@index[3,ac(year_seq),,,,]))); idxCV # a check

# error for observation
Ob.Err <- FLQuant(0,dimnames = list(age=1:9,year=year_seq, iter=1:niter))
Ob.Err[1,,,,,] <- rlnorm(prod(nyear*niter), meanlog = 0, sdlog = 0.05)
Ob.Err[2:9,,,,,] <- Ob.Err[1,,,,,]

### ------------------------------------------------------------------------------------ ###
### --------------------------------------- MSE ---------------------------------------- ###
### ------------------------------------------------------------------------------------ ###

# building catch limit
limit <- data.frame(Year=2015:2020,Item=c(rep("f",6)),Value=c(0.21,0.228,0.092,rep(0.174,3)))


# launch MSE loop

x=1;ii=1
###--------------------------- MSE function ------------------------------ ###
MSE_fun <- function(x){
  # load packages and functions for each cluster
  library(FLCore)
  library(FLash) 
  library(FLAssess) 
  library(FLXSA)
  library(FLBRP)

    # create objects for each iter
  mtf.xsa.obs_nn <- iter(mtf.xsa.obs, x)
  mtf.xsa_nn <- iter(mtf.xsa, x)
  indices.obs_nn <- iterIndicies(indices.obs, x)
  indices_nn <- iterIndicies(indices, x)
  TAC_nn <- iter(TAC, x)
  indices.Err_nn <- iterIndicies(indices.Err, x)
  Ob.Err_nn <- iter(Ob.Err,x)
 
  # iter with year
  for(ii in seq(year_seq)){
    
    # Assess stock from last years catch and indices
    if (ii > 1) {
      mtf.xsa.obs_nn <- assessStock(
        stock = mtf.xsa.obs_nn, index = indices.obs_nn,
        control = xsa.control,
        maxyear = year_seq[ii]-1,
        method="XSA"
      )
    }
    
    # calculate TAC for (year_seq[ii]) using average recruitment 
    ctrl_target <- data.frame(year=(year_seq[ii]),
                              rel.year=NA,
                              quantity=ac(limit$Item[ii]),
                              val=limit$Value[ii],
                              max=NA,min=NA)
    
    ctrl_obj <- fwdControl(ctrl_target)
    stf.xsa <- window(iter(mtf.xsa.obs_nn, 1),start=2014, end=year_seq[ii])
    
    # now fwd the forecast
    stf.xsa <- fwd(stf.xsa, ctrl = ctrl_obj, sr = SR.xsa) 
    
    # record the TAC
    TAC_nn[,ac(year_seq[ii])] <- stf.xsa@catch[,ac(year_seq[ii])]
    
    # advance operational model by one year using TAC (for the real stock)
    mtf.xsa_nn <- advanceStock(
      stock = mtf.xsa_nn, tac = TAC_nn, sr = SR.xsa,
      maxyear = year_seq[ii]-1,
      sr.residuals = iter(trim(SR.Err, year=year_seq[ii]), x),
      sr.residuals.mult = TRUE
    )      
    
    # record observed catches
    mtf.xsa.obs_nn <- observeCatch(stockReal=mtf.xsa_nn, stockObs=mtf.xsa.obs_nn, obsErr = Ob.Err_nn,
                                   minyear=NULL, maxyear=NULL)
    
    # record real and observed indices
    indices_nn <- observeIndex(stock=mtf.xsa_nn, index=indices_nn, q.hat=q.hat, obsErr=NULL, minyear = year_seq[1])
    indices.obs_nn <- observeIndex(stock=mtf.xsa_nn, index=indices.obs_nn, q.hat=q.hat, obsErr=indices.Err, 
                                   minyear = year_seq[1])
    
    # add penalty when SSB is between 7500 and 12000
    if (ssb(mtf.xsa_nn)[,ac(year_seq[ii]),,,,]>7500 &
        ssb(mtf.xsa_nn)[,ac(year_seq[ii]),,,,]<12500 )  {
      stock.n(mtf.xsa_nn)[,ac(year_seq[ii]),,,,] <- stock.n(mtf.xsa_nn)[,ac(year_seq[ii]),,,,]/1.027757
    }
  }
  
  mtf.xsa.obs_nn <- assessStock(
    stock = mtf.xsa.obs_nn, index = indices.obs_nn,
    control = xsa.control,
    maxyear = year_seq[ii],
    method="XSA"
  )
  
  # recorded the finished iteration
  RES <- list(
    mtf.xsa.obs_nn = mtf.xsa.obs_nn,
    mtf.xsa_nn = mtf.xsa_nn
  )
  
  sink(file="output.txt", append = TRUE)
  print(paste("iter", x, "completed @", Sys.time()))
  sink()
  
  return(RES)
}  # end MSE


###----------------- run MSE with cluster and export the result --------------- ###
cl <- makeCluster(ncores, type="SOCK", outfile="tmp.txt")
clusterExport(cl, list=ls(), envir = .GlobalEnv)

t1 <- Sys.time()
res <- parLapply(cl, 1:niter, function(x) try(MSE_fun(x), silent= T) )
stopCluster(cl)
Sys.time()- t1

# Save the results
save(res, file="results/Projection with F Gen.SR scaled.Rdata")




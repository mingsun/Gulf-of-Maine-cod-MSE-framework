

library(FLCore)
library(FLash)
library(FLBRP)
library(FLAssess)
library(FLXSA)
library(snow)

load(file = "R/data/survey.hind.R")  # indices
load(file = "R/data/gom.hind.R")     # the hindsight gom file
load(file = "R/data/q.hat.R")        # q.hat
load(file = "R/data/SRR.1st.R")      # SR
load(file = "R/data/SRR.2nd.R")      # SR
load(file = "R/data/SRR.3rd.R")      # SR

# true data
load(file = "R/data/gom.xsa.R")
load(file = "R/data/survey.R")

source("R/functions/functions for GoM cod.R")
gom <- gom.hind
indices <- indices.hind

##############################################################################################
#### --------------------------------- General Setting ---------------------------------- ####
##############################################################################################

###--------------------------- Fundemental Parameters ------------------------------ ###

niter <- 200
ncores <- 1 #as.numeric(Sys.getenv("NUMBER_OF_PROCESSORS"))
nn <- split(1:niter, 1:niter)

year_seq <- seq(range(gom)["minyear"],range(gom)["maxyear"]) # starting from 1982
nyear <- length(year_seq)

###--------------------------- Input data Processing ------------------------------ ###

mtf.xsa <- gom
indices <- window(indices, end=tail(year_seq,1))

# expand to niter
mtf.xsa <- propagate(mtf.xsa, niter)
for(i in seq(indices)){
  indices[[i]] <- propagate(indices[[i]], niter)
}


TAC <- catch(mtf.xsa)
catch.input <- read.csv("R/data/catch.input.csv",row.names = 1)


###--------------------------- Stochasticity ------------------------------ ###

set.seed(1)
# error for flux
Flux.Err <- FLQuant(rnorm(niter,mean = 1.013205,sd= 0.5),
                    dimnames = list(year=year_seq, iter=1:niter))


##############################################################################################
#### --------------------------------------- MSE ---------------------------------------- ####
##############################################################################################

# building scenarios for different levels of flux variation
for (ref in seq(0.01,0.04,0.01)) {
  
  # make duplicates of objects
  mtf.xsa_dup     <- mtf.xsa
  
  Flux.Err <- FLQuant(rnorm(niter,mean = 1.013205,sd= ref),
                      dimnames = list(year=year_seq, iter=1:niter))
  
  ###--------------------------- MSE function ------------------------------ ###
  MSE_fun <- function(x){
    library(FLCore)
    library(FLash) 
    library(FLAssess) 
    library(FLXSA)
    library(FLBRP)
    
    # create objects for each iter
    mtf.xsa_nn <- iter(mtf.xsa, x)
    TAC_nn <- iter(TAC, x)
    
    for(ii in 2:(length(year_seq))){ # starting from 1982
      
      # update catch time series
      TAC_nn[,ac(year_seq[ii]),,,,] <- as.numeric(catch.input$Catch[ii])
      
      if (year_seq[ii]<=1988){SR <- SR.1st}
      if (year_seq[ii]>1988 & year_seq[ii]<=2002){SR <- SR.2nd}
      if (year_seq[ii]>2002){SR <- SR.3rd}

      # advance operational model by one year using TAC (for the real stock)
      mtf.xsa_nn <- advanceStock(
        stock = mtf.xsa_nn, tac = TAC_nn, sr = SR,
        maxyear = year_seq[ii]-1, 
        sr.residuals = exp(SR@residuals)[,ac(year_seq[ii]),,,,],
        sr.residuals.mult = TRUE
      )
      
      # add penalty when SSB is between 7500 and 12000
      if (ssb(mtf.xsa_nn)[,ac(year_seq[ii]),,,,]>7500 &
          ssb(mtf.xsa_nn)[,ac(year_seq[ii]),,,,]<12500)  {
        stock.n(mtf.xsa_nn)[,ac(year_seq[ii]),,,,] <- 
          stock.n(mtf.xsa_nn)[,ac(year_seq[ii]),,,,]/as.numeric(Flux.Err[,ac(year_seq[ii]),,,,x])
      }
      
    } 
    
    sink(file="output.txt", append = TRUE)
    print(paste("Scenario",ref, "iter", x, "completed @", Sys.time()))
    sink()
    
    return(mtf.xsa_nn)
  } # end MSE
  
  ###----------------- run MSE with cluster and export the result --------------- ###
  cl <- makeCluster(ncores, type="SOCK", outfile="tmp.txt")
  clusterExport(cl, list=ls(), envir = .GlobalEnv)
  
  t1 <- Sys.time()
  res <- parLapply(cl, 1:niter, function(x) try(MSE_fun(x), silent= T) )
  stopCluster(cl)
  Sys.time()- t1
  
  save( res,   file=paste("results/Hindsight MSE Seg.SR/Ref ",  ref,".Rdata", sep="") )
  
  # restore objects
  mtf.xsa     <- mtf.xsa_dup
}


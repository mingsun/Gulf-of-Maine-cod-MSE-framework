

library(FLCore)
library(FLash)
library(FLBRP)
library(FLAssess)
library(ggplot2)

load(file = "R/data/survey.hind.R")  # indices
load(file = "R/data/gom.hind.R")     # the hindsight gom file
load(file = "R/data/SRR.xsa.R")      # SR
load(file = "R/data/q.hat.R")        # q.hat

# true data
load(file = "R/data/gom.xsa.R")
load(file = "R/data/survey.R")

source("R/functions/functions for GoM cod.R")
result_path <- "results/validation M 0.2/"
gom <- gom.hind
indices <- indices.hind

### ------------------------------------------------------------------------------------ ###
### --------------------------------- General Setting ---------------------------------- ###
### ------------------------------------------------------------------------------------ ###

###--------------------------- Fundemental Parameters ------------------------------ ###

year_seq <- seq(range(gom)["minyear"],range(gom)["maxyear"]) # starting from 1982
nyear <- length(year_seq)

###--------------------------- Input data Processing ------------------------------ ###

mtf.xsa <- gom
indices <- window(indices, end=tail(year_seq,1))

TAC <- catch(mtf.xsa)
catch.input <- read.csv("R/data/catch.input.csv",row.names = 1)

### ------------------------------------------------------------------------------------ ###
### --------------------------------------- MSE ---------------------------------------- ###
### ------------------------------------------------------------------------------------ ###


### ------------- 3. with one SR and penalty ------------ ###

# create test objects
mtf.xsa_nn <- mtf.xsa
indices_nn <- indices
TAC_nn <- TAC


# optimize the minimum residual with flux penalty
min.RSS <- function(par) {
  
  for(ii in 2:(length(year_seq))){ # starting from 1982
    # update catch time series
    TAC_nn[,ac(year_seq[ii]),,,,] <- as.numeric(catch.input$Catch[ii])
    SR <- SR.xsa
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
      stock.n(mtf.xsa_nn)[,ac(year_seq[ii]),,,,] <- stock.n(mtf.xsa_nn)[,ac(year_seq[ii]),,,,]/par
    }
  } 
  
  sum((ssb(mtf.xsa_nn)-ssb(gom.xsa))^2)/33 # value to be minimised
}

result <- optim(par=1,fn = min.RSS)

### ------------- 4. with different SRRs ------------ ###

load(file = "R/data/SRR.1st.R")      # SR
load(file = "R/data/SRR.2nd.R")      # SR
load(file = "R/data/SRR.3rd.R")      # SR

# create test objects
mtf.xsa_nn <- mtf.xsa
indices_nn <- indices
TAC_nn <- TAC


# optimize the minimum residual with flux penalty
min.RSS <- function(par) {
  
  for(ii in 2:(length(year_seq))){ # starting from 1982
    # update catch time series
    TAC_nn[,ac(year_seq[ii]),,,,] <- as.numeric(catch.input$Catch[ii])
    
    if (year_seq[ii]<=1988){SR <- SR.1st}
    if (year_seq[ii]>1988 & year_seq[ii]<=2002){SR <- SR.2nd}
    # if (year_seq[ii]<=2002){SR <- SR.xsa}
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
      stock.n(mtf.xsa_nn)[,ac(year_seq[ii]),,,,] <- stock.n(mtf.xsa_nn)[,ac(year_seq[ii]),,,,]/par
    }
  } 
  
  sum((ssb(mtf.xsa_nn)-ssb(gom.xsa))^2)/33 # value to be minimised
}

result.2 <- optim(par=1, fn = min.RSS)





library(FLCore)
library(FLash)
library(FLBRP)
library(FLAssess)
library(ggplot2)

load(file = "R/data/survey.hind.R")  # indices
load(file = "R/data/gom.hind.R")     # the hindsight gom file
load(file = "R/data/SRR.1st.R")      # SR
load(file = "R/data/SRR.2nd.R")      # SR
load(file = "R/data/SRR.3rd.R")      # SR
load(file = "R/data/q.hat.R")        # q.hat

# true data
load(file = "R/data/gom.xsa.R")
load(file = "R/data/survey.R")

source("R/functions/functions for GoM cod.R")
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

### ------------- 1. with no error ------------ ###

# create test objects
mtf.xsa_nn <- mtf.xsa
indices_nn <- indices
TAC_nn <- TAC

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
    sr.residuals.mult = FALSE
  )
} 

no.error.seg.SR <- ssb(mtf.xsa_nn)
save(no.error.seg.SR,file="results/Deterministic MSE SSB/no error seg.SR.Rdata")


### ------------- 2. with actual error ------------ ###

# create test objects
mtf.xsa_nn <- mtf.xsa
indices_nn <- indices
TAC_nn <- TAC

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
} 

ggplot() +
  geom_line(aes(x=1982:2014,y=as.numeric(ssb(mtf.xsa_nn))),color="black",linetype=2,size=1) +
  geom_line(aes(x=1982:2014,y=as.numeric(ssb(gom.xsa))),color="red",size=1)

with.error.seg.SR<- ssb(mtf.xsa_nn)
save(with.error.seg.SR,file="results/Deterministic MSE SSB/with error seg.SR.Rdata")


### ------------- 3. with actual error and SSB penalty ------------ ###

# create test objects
mtf.xsa_nn <- mtf.xsa
indices_nn <- indices
TAC_nn <- TAC

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
      ssb(mtf.xsa_nn)[,ac(year_seq[ii]),,,,]<12500 )  {
    stock.n(mtf.xsa_nn)[,ac(year_seq[ii]),,,,] <- stock.n(mtf.xsa_nn)[,ac(year_seq[ii]),,,,]/1.013205
    # stock.n(mtf.xsa_nn)[,ac(year_seq[ii]),,,,] <- stock.n(mtf.xsa_nn)[,ac(year_seq[ii]),,,,]/1.027757
  }
  
  # record real and observed indices
  indices_nn <- observeIndex(stock=mtf.xsa_nn, index=indices_nn, q.hat=q.hat, obsErr=NULL, minyear = year_seq[1])
} 

ggplot() +
  geom_line(aes(x=1982:2014,y=as.numeric(ssb(mtf.xsa_nn))),color="black",linetype=2,size=1) +
  geom_line(aes(x=1982:2014,y=as.numeric(ssb(gom.xsa))),color="red",size=1)

ssb.with.seg.SR <- ssb(mtf.xsa_nn)
save(ssb.with.seg.SR,file="results/ssb.with.seg.SR.Rdata")

indices.with.seg.SR <- indices_nn
save(indices.with.seg.SR,file="results/indices.with.seg.SR.Rdata")

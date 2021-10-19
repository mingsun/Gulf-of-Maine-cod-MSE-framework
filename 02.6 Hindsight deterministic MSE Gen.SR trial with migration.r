

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

mig <- read.csv("R/data/Migration.csv")[,1:3]
stock.n(gom)

### ------------------------------------------------------------------------------------ ###
### --------------------------------------- MSE ---------------------------------------- ###
### ------------------------------------------------------------------------------------ ###

### ------------- 1. with no error ------------ ###

# create test objects
mtf.xsa_nn <- mtf.xsa
indices_nn <- indices
TAC_nn <- TAC

ii=2
ii=ii+1

for(ii in 2:(length(year_seq))){ # starting from 1982
  
  # update catch time series
  TAC_nn[,ac(year_seq[ii]),,,,] <- as.numeric(catch.input$Catch[ii])
  SR <- SR.xsa
  
  # remove the migration influence
  migrated <- mtf.xsa_nn
  migrated.scale <- as.numeric((apply(stock.n(migrated)[,ac(year_seq[ii-1]),,,,],2,sum)/1000 + mig$migrants.in.million[ii-1])/
                                 (apply(stock.n(migrated)[,ac(year_seq[ii-1]),,,,],2,sum)/1000))
  stock.n(migrated)[,ac(year_seq[ii-1]),,,,] <- stock.n(migrated)[,ac(year_seq[ii-1]),,,,]*migrated.scale
  
  # advance with catch for age 2:9
  migrated <- advanceStock(
    stock = migrated, tac = TAC_nn, sr = SR,
    maxyear = year_seq[ii]-1, 
    sr.residuals = exp(SR@residuals)[,ac(year_seq[ii]),,,,], 
    sr.residuals.mult = FALSE
  )
  
  # advance for recruitment
  mtf.xsa_nn <- advanceStock(
    stock = mtf.xsa_nn, tac = TAC_nn, sr = SR,
    maxyear = year_seq[ii]-1, 
    sr.residuals = exp(SR@residuals)[,ac(year_seq[ii]),,,,], 
    sr.residuals.mult = FALSE
  )
  
  stock.n(mtf.xsa_nn)[2:9,ac(year_seq[ii]),,,,] <- stock.n(migrated)[2:9,ac(year_seq[ii]),,,,]/migrated.scale
  
} 

ggplot() +
  geom_line(aes(x=1982:2014,y=as.numeric(ssb(mtf.xsa_nn))),color="black",linetype=2,size=1) +
  geom_line(aes(x=1982:2014,y=as.numeric(ssb(gom.xsa))),color="red",size=1)

ssb.without.error<- ssb(mtf.xsa_nn)
save(ssb.without.error,file="results/ssb.without.error.Rdata")


### ------------- 2. with actual error ------------ ###

# create test objects
mtf.xsa_nn <- mtf.xsa
indices_nn <- indices
TAC_nn <- TAC

for(ii in 2:(length(year_seq))){ # starting from 1982
  
  # update catch time series
  TAC_nn[,ac(year_seq[ii]),,,,] <- as.numeric(catch.input$Catch[ii])
  SR <- SR.xsa
  
  # remove the migration influence
  migrated <- mtf.xsa_nn
  migrated.scale <- as.numeric((apply(stock.n(migrated)[,ac(year_seq[ii-1]),,,,],2,sum)/1000 + mig$migrants.in.million[ii-1])/
                                 (apply(stock.n(migrated)[,ac(year_seq[ii-1]),,,,],2,sum)/1000))
  stock.n(migrated)[,ac(year_seq[ii-1]),,,,] <- stock.n(migrated)[,ac(year_seq[ii-1]),,,,]*migrated.scale
  
  # advance with catch for age 2:9
  migrated <- advanceStock(
    stock = migrated, tac = TAC_nn, sr = SR,
    maxyear = year_seq[ii]-1, 
    sr.residuals = exp(SR@residuals)[,ac(year_seq[ii]),,,,], 
    sr.residuals.mult = TRUE
  )
  
  # advance for recruitment
  mtf.xsa_nn <- advanceStock(
    stock = mtf.xsa_nn, tac = TAC_nn, sr = SR,
    maxyear = year_seq[ii]-1, 
    sr.residuals = exp(SR@residuals)[,ac(year_seq[ii]),,,,], 
    sr.residuals.mult = TRUE
  )
  
  stock.n(mtf.xsa_nn)[2:9,ac(year_seq[ii]),,,,] <- stock.n(migrated)[2:9,ac(year_seq[ii]),,,,]
} 

ggplot() +
  geom_line(aes(x=1982:2014,y=as.numeric(ssb(mtf.xsa_nn))),color="black",linetype=2,size=1) +
  geom_line(aes(x=1982:2014,y=as.numeric(ssb(gom.xsa))),color="red",size=1)

ssb.with.error<- ssb(mtf.xsa_nn)
save(ssb.with.error,file="results/ssb.with.error.Rdata")


### ------------- 3. with actual error and SSB penalty ------------ ###

# create test objects
mtf.xsa_nn <- mtf.xsa
indices_nn <- indices
TAC_nn <- TAC

for(ii in 2:(length(year_seq))){ # starting from 1982
  
  # update catch time series
  TAC_nn[,ac(year_seq[ii]),,,,] <- as.numeric(catch.input$Catch[ii])
  SR <- SR.xsa
  
  # remove the migration influence
  migrated <- mtf.xsa_nn
  migrated.scale <- as.numeric((apply(stock.n(migrated)[,ac(year_seq[ii-1]),,,,],2,sum)/1000 + mig$migrants.in.million[ii-1])/
                                 (apply(stock.n(migrated)[,ac(year_seq[ii-1]),,,,],2,sum)/1000))
  stock.n(migrated)[,ac(year_seq[ii-1]),,,,] <- stock.n(migrated)[,ac(year_seq[ii-1]),,,,]*migrated.scale
  
  # advance with catch for age 2:9
  migrated <- advanceStock(
    stock = migrated, tac = TAC_nn, sr = SR,
    maxyear = year_seq[ii]-1, 
    sr.residuals = exp(SR@residuals)[,ac(year_seq[ii]),,,,], 
    sr.residuals.mult = TRUE
  )
  
  # advance for recruitment
  mtf.xsa_nn <- advanceStock(
    stock = mtf.xsa_nn, tac = TAC_nn, sr = SR,
    maxyear = year_seq[ii]-1, 
    sr.residuals = exp(SR@residuals)[,ac(year_seq[ii]),,,,], 
    sr.residuals.mult = TRUE
  )
  
  stock.n(mtf.xsa_nn)[1,ac(year_seq[ii]),,,,] <- stock.n(migrated)[1,ac(year_seq[ii]),,,,]
  
  # add penalty when SSB is between 7500 and 12000
  if (ssb(mtf.xsa_nn)[,ac(year_seq[ii]),,,,]>7500 &
      ssb(mtf.xsa_nn)[,ac(year_seq[ii]),,,,]<12500)  {
    stock.n(mtf.xsa_nn)[,ac(year_seq[ii]),,,,] <- stock.n(mtf.xsa_nn)[,ac(year_seq[ii]),,,,]/1.027757
  }
} 

ggplot() +
  geom_line(aes(x=1982:2014,y=as.numeric(ssb(mtf.xsa_nn))),color="black",linetype=2,size=1) +
  geom_line(aes(x=1982:2014,y=as.numeric(ssb(gom.xsa))),color="red",size=1)

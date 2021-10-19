
library(reshape2)
library(FLCore)
   


### ------------------------------------------------------------------------------------ ###
### ------------------------------ Erase Known Stock Info ------------------------------ ###
### ------------------------------------------------------------------------------------ ###

load(file = "R/data/gom.xsa.R") ; gom.hind <- gom.xsa
erase.year <- 1983:2014

gom.hind@catch[,ac(erase.year),,,,] <- gom.hind@catch.n[,ac(erase.year),,,,] <- gom.hind@landings[,ac(erase.year),,,,] <- 
  gom.hind@landings.n[,ac(erase.year),,,,] <- gom.hind@stock.n[,ac(erase.year),,,,] <- gom.hind@harvest[,ac(erase.year),,,,] <-  NA

gom.hind@harvest[,ac(erase.year),,,,] <- catch.sel(gom.xsa)[,ac(erase.year),,,,] 

save(gom.hind,file="R/data/gom.hind.R")

### ------------------------------------------------------------------------------------ ###
### ----------------------------- Erase Known Survey Infor ----------------------------- ###
### ------------------------------------------------------------------------------------ ###

load(file="R/data/survey.R"); indices.hind <- indices
erase.year <- 1983:2014

indices.hind[[1]]@index[,ac(erase.year),,,,] <- indices.hind[[1]]@catch.n[,ac(erase.year),,,,] <- 
  indices.hind[[2]]@index[,ac(erase.year),,,,] <- indices.hind[[2]]@catch.n[,ac(erase.year),,,,] <-
     indices.hind[[3]]@index[,ac(erase.year),,,,] <- indices.hind[[3]]@catch.n[,ac(erase.year),,,,] <- NA

save(indices.hind,file="R/data/survey.hind.R")

### ------------------------------------------------------------------------------------ ###
### ------------------------------ Extract Observed Data ------------------------------- ###
### ------------------------------------------------------------------------------------ ###

load(file = "R/data/gom.xsa.R")   

# Total Catch Data 
Catch <- as.numeric(gom.xsa@catch)

# Catch number at age
CAA <- as.data.frame(gom.xsa@catch.n)
CAA <- dcast(CAA,year ~ age, value.var = "data")
colnames(CAA)[1] <- c("Year")
colnames(CAA)[2:10] <- c(paste("Age",colnames(CAA)[2:10]))

CAA$Catch <- Catch
catch.input <- CAA;remove(CAA)

write.csv(catch.input,"R/data/catch.input.csv")

### ------------------------------------------------------------------------------------ ###
### -------------------------------------- M-model ------------------------------------- ###
### ------------------------------------------------------------------------------------ ###

# M=0.2 model all set
# M 0.2-0.4 model; 0.2 (1982-1988),0.4 (2003-2011), with a linear ramp in between

a <- 0.2/(2003-1988)
a*1988+b=0.2
b <- 0.2-a*1988

M.input <- data.frame(Year=1982:2014,M=round(c(rep(0.2,7),a*c(1989:2002)+b,rep(0.4,12)),3))
  
write.csv(M.input,"R/data/M.input.csv")                    

### ------------------------------------------------------------------------------------ ###
### ----------------------- Pre-define Stochasticity distribution ---------------------- ###
### ------------------------------------------------------------------------------------ ###



# recruitment stochasticity Monte Carlo simulation (ask Lisa)


# sensitivity analysis with bias in catch and indices observation, bias in M variation, bias in implementation (with both F and Catch)




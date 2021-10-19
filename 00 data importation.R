
setwd("C://Users/sunming1991/Desktop/phd/paper/data-limited vs data-rich/")

library(FLCore)

### ---------------------------------- load data ---------------------------------- ###
#####################   stock

catch.n <- read.csv("data from Chen/reorganize GOM data/catch data 1.csv",sep="",header=FALSE)
catch.n.matrix <- t(as.matrix(catch.n[,1:9]))
catch.n.flq <- FLQuant(catch.n.matrix, dimnames=list(age=1:9, year = 1982:2014),quant="age")

catch.matrix <- t(as.matrix(catch.n[,10]))
catch.flq <- FLQuant(catch.matrix, dimnames=list(year = 1982:2014),quant="age")

weight.age1 <- read.csv("data from Chen/reorganize GOM data/weight matrix 1.csv",sep="",header=FALSE)  # for catch calculation
weight.age1.matrix <- t(as.matrix(weight.age1))
weight.age1.flq <- FLQuant(weight.age1.matrix, dimnames=list(age=1:9, year = 1982:2014),quant="age")

weight.age2 <- read.csv("data from Chen/reorganize GOM data/weight matrix 2.csv",sep="",header=FALSE)  # for biomass calculation
weight.age2.matrix <- t(as.matrix(weight.age2))
weight.age2.flq <- FLQuant(weight.age2.matrix, dimnames=list(age=1:9, year = 1982:2014),quant="age")

discard <- read.csv("data from Chen/reorganize GOM data/discard.csv",header=FALSE)
discard.matrix <- t(as.matrix(discard[,3]))
discard.flq <- FLQuant(discard.matrix, dimnames=list(year = 1982:2014),quant="age")

maturity <- read.csv("data from Chen/reorganize GOM data/Maturity.csv",sep="",header=FALSE)
maturity.matrix <- t(as.matrix(maturity))
maturity.flq <- FLQuant(maturity.matrix, dimnames=list(age=1:9, year = 1982:2014),quant="age")


apply(catch.n[,1:9]*weight.age1,1,sum)

# [1] "catch"        "catch.n"      "catch.wt"     "discards"     "discards.n"   "discards.wt"  "landings"     "landings.n"   "landings.wt" 
# [10] "stock"        "stock.n"      "stock.wt"     "m"            "mat"          "harvest"      "harvest.spwn" "m.spwn"       "name"        
# [19] "desc"         "range" 
data(ple4)

cod.fishery <- ple4
cod.fishery@catch <-   catch.flq
cod.fishery@catch.n <- catch.n.flq
cod.fishery@catch.wt <- weight.age1.flq
cod.fishery@discards <- FLQuant(0,dimnames=list(year = 1982:2014),quant="age")
cod.fishery@discards.n <- FLQuant(0,dimnames=list(age=1:9, year = 1982:2014),quant="age")
cod.fishery@discards.wt <- cod.fishery@catch.wt
cod.fishery@landings <- cod.fishery@catch
cod.fishery@landings.n <- cod.fishery@catch.n
cod.fishery@landings.wt <- cod.fishery@catch.wt
cod.fishery@stock <- FLQuant(0,dimnames=list(age=1:9, year = 1982:2014),quant="age")
cod.fishery@stock.n <- FLQuant(0,dimnames=list(age=1:9, year = 1982:2014),quant="age")
cod.fishery@stock.wt <- weight.age2.flq
cod.fishery@m <- FLQuant(0.2,dimnames=list(age=1:9, year = 1982:2014),quant="age")
cod.fishery@mat <- maturity.flq
cod.fishery@harvest <- FLQuant(0,dimnames=list(age=1:9, year = 1982:2014),quant="age")
cod.fishery@harvest.spwn <- FLQuant(0.25,dimnames=list(age=1:9, year = 1982:2014),quant="age")
cod.fishery@m.spwn <- FLQuant(0.25,dimnames=list(age=1:9, year = 1982:2014),quant="age")
cod.fishery@name <- "Gulf of Maine cod"
cod.fishery@desc <- "NA"
cod.fishery@range  
range(cod.fishery)[c("min","max","plusgroup","minyear","maxyear","minfbar","maxfbar")] <- as.numeric(c(1,9,9,1982,2014,5,5))

units(catch(cod.fishery)) <- units(discards(cod.fishery)) <- units(landings(cod.fishery)) <- units(stock(cod.fishery)) <- 'tonnes'
units(catch.n(cod.fishery)) <- units(discards.n(cod.fishery)) <- units(landings.n(cod.fishery)) <- units(stock.n(cod.fishery)) <- '1000'
units(catch.wt(cod.fishery)) <- units(discards.wt(cod.fishery)) <- units(landings.wt(cod.fishery)) <- units(stock.wt(cod.fishery)) <- 'kg'
units(harvest(cod.fishery)) <- 'f'

save(cod.fishery,file="R/data/GoMcod.R")

#####################   survey indices
index.1 <- read.csv("data from Chen/reorganize GOM data/index data 1.csv",header=FALSE)
index.1.matrix <- t(as.matrix(index.1[,4:12]))
index.1.flq <- FLQuant(index.1.matrix, dimnames=list(age=1:9, year = 1982:2014))
indices.1 <- FLIndex(index = index.1.flq,catch.n=index.1.flq,effort=FLQuant(1, dimnames=list(year = 1982:2014)),
                    name="NEFSCspring")
indices.1@sel.pattern <- FLQuant(matrix(c(0.05,0.2,0.4,0.79,0.9,1,1,1,1),nrow = 9,ncol = 33),
                                 dimnames=list(age=1:9, year = 1982:2014))
range(indices.1)[c('startf', 'endf')] <- c(3/12,5/12)


index.2 <- read.csv("data from Chen/reorganize GOM data/index data 2.csv",header=FALSE)
index.2.matrix <- t(as.matrix(index.2[,4:12]))
index.2.flq <- FLQuant(index.2.matrix, dimnames=list(age=1:9, year = 1982:2014))
indices.2 <- FLIndex(index = index.2.flq,catch.n=index.2.flq,effort=FLQuant(1, dimnames=list(year = 1982:2014)),
                     name="NEFSCfall")
indices.2@sel.pattern <- FLQuant(matrix(c(0.05,0.2,0.4,0.79,0.9,1,1,1,1),nrow = 9,ncol = 33),
                                 dimnames=list(age=1:9, year = 1982:2014))
range(indices.2)[c('startf', 'endf')] <- c(9/12,11/12)


index.3 <- read.csv("data from Chen/reorganize GOM data/index data 3.csv",header=FALSE)
index.3.matrix <- t(as.matrix(index.3[,4:12]))
index.3.flq <- FLQuant(index.3.matrix, dimnames=list(age=1:9, year = 1982:2014))
indices.3 <- FLIndex(index = index.3.flq,catch.n=index.3.flq,effort=FLQuant(1, dimnames=list(year = 1982:2014)),
                     name="MAspring")
indices.3@sel.pattern <- FLQuant(matrix(c(1,0.8,0.6,0.5,0.4,0.2,0,0,0),nrow = 9,ncol = 33),
                                 dimnames=list(age=1:9, year = 1982:2014))
range(indices.3)[c('startf', 'endf')] <- c(3/12,5/12)

indices <- FLIndices(indices.1,indices.2,indices.3)
save(indices,file="R/data/survey.R")




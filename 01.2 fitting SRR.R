

library(FLCore)
library(FLAssess)
library(FLXSA)
library(FLBRP)

load(file = "R/data/gom.xsa.R")

##############################################################################################
#### ----------------------------------- Fitting SRR ------------------------------------ ####
##############################################################################################

# SRR with xsa
# gom.new <- trim(gom.xsa, year=c(1982:1987,1989:2009)) #recruitment in 1988 is extremly high
gom.new <- trim(gom.xsa, year=c(1982:2014))
SR.xsa <- as.FLSR(gom.new)
model(SR.xsa) <- ricker()
SR.xsa<-fmle(SR.xsa,fixed = list(a = 0.7,b=2.98*10^(-5)))
params(SR.xsa)
# SR.xsa@params[[1]] <- 0.7
# SR.xsa@params[[2]] <- 2.98*10^(-5)
plot(SR.xsa)

save(SR.xsa,file="R/data/SRR.xsa.R")

# -------- 1982:1989
gom <- gom.xsa

gom.new <- trim(gom, year=c(1982:1988))           # 4.72078290 0.00010604
SR.1st <- as.FLSR(gom.new)
model(SR.1st) <- ricker
SR.1st<-fmle(SR.1st)
params(SR.1st)

plot(SR.1st)
save(SR.1st,file="R/data/SRR.1st.R")


# --------- 1989:2002
gom.new <- trim(gom, year=c(1988:2002))           # 0.60752 2.9033e-05
SR.2nd <- as.FLSR(gom.new)
model(SR.2nd) <- ricker
SR.2nd<-fmle(SR.2nd)
params(SR.2nd)

plot(SR.2nd)
save(SR.2nd,file="R/data/SRR.2nd.R")

# --------- 2003:2014
gom.new <- trim(gom, year=c(2002:2014))           # 0.9439108 0.0001318
SR.3rd <- as.FLSR(gom.new)
model(SR.3rd) <- ricker
SR.3rd<-fmle(SR.3rd)
params(SR.3rd)

plot(SR.3rd)
save(SR.3rd,file="R/data/SRR.3rd.R")


# --------- 1989:2014
gom.new <- trim(gom, year=c(1988:2014))           # 0.60752 2.9033e-05
SR.23 <- as.FLSR(gom.new)
model(SR.23) <- ricker
SR.23<-fmle(SR.23)
params(SR.23)

plot(SR.23)
save(SR.23,file="R/data/SRR.23.R")


##############################################################################################
#### ---------------------------- Initial Reference points ------------------------------ ####
##############################################################################################




temp.prp <- brp(FLBRP(gom.new, SR.xsa))
dimnames(refpts(temp.prp))$refpt[6] <- c('spr.40')
temp.prp <- brp(FLBRP(temp.prp))
temp.refp <- refpts(temp.prp)
temp.refp

refp <- as.data.frame(temp.refp)

write.csv(refp,file="R/data/collapsing rfp.csv")


##############################################################################################
#### ---------------------------- sel at age ------------------------------ ####
##############################################################################################

sel.age <- c(0.006584904,0.074039558,0.318910453,0.769927456,1.000000000,
             0.992274604,0.960866749,0.763601010,0.763601010)
save(sel.age,file="R/data/sel.age.R")

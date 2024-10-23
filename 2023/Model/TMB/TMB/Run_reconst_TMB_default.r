#'==============================================================================
#   Kuskokwim Chinook Salmon Run Reconstruction R-code TMB Implementation 
#   Written by  H. Hamaazaki 
#   Date  06/14/2023         
#   This reads TMB.ccp 
#'==============================================================================

#'==============================================================================
#  1.0  Initialize working Environment ----                                        
#'==============================================================================
library(TMB)
#library(TMBdebug)
rm(list=ls(all=TRUE))
# Set Base directry 
Base_dir <- file.path('C:','Projects','Kuskokwim_River','Chinook_reconst','Hamachan','TMB')
# Specify TMB directry 
TMBdir <- file.path(Base_dir,'Model')
# Specify Data directory 
Data_dir <- file.path(Base_dir,'Data')
# Install ADMB read write R source code:
source(file.path('C:','Projects','Statistics_Worksheets','R_functions','ADMBtoR.r'))


# Specify TMB Source file name 
TMB.model <- 'Chinook_reconst_TMB_default.cpp'
# Spedify Input CSV data file

# Specify the name of accompanying program file 
datafn <- file.path(Data_dir,'Kusko_Chinook_2023.dat')
#  Read ADMB data  to list 
TMB.data <- datatoR(datafn) 
TMB.data$w_esc <- t(TMB.data$esc_w)
TMB.data$a_esc <- t(TMB.data$esc_a)


cpue <- with(TMB.data ,ifelse(ceff>0,ccat/ceff,0))
# Unrestricted Catch 
ureg <- with(TMB.data ,ifelse(creg==1,1,0))
ur <- rowSums(ureg*cpue,na.rm=TRUE)
urp <- rowSums(ureg*TMB.data$testf,na.rm=TRUE)

# Restricted Catch 
ureg <- with(TMB.data ,ifelse(creg==2,1,0))
rst <- rowSums(ureg*cpue,na.rm=TRUE)
rp <- rowSums(ureg*TMB.data $testf,na.rm=TRUE)

# Monofilament Catch 
ureg <- with(TMB.data ,ifelse(creg==3|creg==5,1,0))
nst <- rowSums(ureg*cpue,na.rm=TRUE)
np <- rowSums(ureg*TMB.data$testf,na.rm=TRUE)

TMB.data$testp <-  as.matrix(t(cbind(urp,rp,np)))
TMB.data$cpue <-  as.matrix(t(cbind(ur,rst,nst)))

parameters <- list(log_trun=rep(12.5,TMB.data$nyear), log_wesc=rep(5.0,TMB.data$nweir), log_aesc = rep(4.0,TMB.data$naerial), 
log_cvw = 0.5, log_cva=0.5,log_q = rep(-11,2), log_cvq = 0.5)

p.lower <- c(log_trun=TMB.data$minrun, log_wesc=rep(0.0,TMB.data$nweir), log_aesc = rep(0.0,TMB.data$naerial), 
log_cvw = -10.0, log_cva=-10.0,log_q = rep(-12,2), log_cvq = -10.0)

p.upper <- c(log_trun=TMB.data$ubrun, log_wesc=rep(7.0,TMB.data$nweir), log_aesc = rep(7.0,TMB.data$naerial),log_cvw = 1.0, log_cva=1.0, 
log_q = rep(-9,2), log_cvq = 1.0)
#'###########################################################################
# Run the ADMB ----
#'###########################################################################	
# Compile model
setwd(TMBdir)
compile(TMB.model,flags="-O0 -g",DLLFLAGS="",libtmb=FALSE)

dyn.load(dynlib("Chinook_reconst_TMB_default"))
model <- MakeADFun(TMB.data, parameters, DLL="Chinook_reconst_TMB_default", silent=FALSE) 
fit <- nlminb(model$par, model$fn, model$gr, control = list(eval.max = 10000, iter.max=10000), upper=p.upper, lower=p.lower)
rep <- sdreport(model)
print(summary(rep))

par1 <-summary(rep)


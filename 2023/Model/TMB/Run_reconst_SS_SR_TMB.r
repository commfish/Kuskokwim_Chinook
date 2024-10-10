#'==============================================================================
#   Kuskokwim Chinook Salmon Run Reconstruction R-code TMB Implementation 
#
#   Written by  H. Hamaazaki 
#   Date  06/14/2023         
#   
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
TMB.model <- 'Chinook_reconst_SS_SR_TMB.cpp'
# Spedify Input CSV data file
data_file1 = 'Kusko_RR_Input_Oct_2022.csv'

kusko.data <- read.csv(file.path(Data_dir,data_file1),header=TRUE, na.string='')

# Specify the name of accompanying program file 
datafn <- file.path(Data_dir,'Kusko_Chinook_22_new_v.dat')
#  Read ADMB data  to list 
TMB.data <- datatoR(datafn) 


cpue <- with(TMB.data ,ccat/ceff)
# Unrestricted Catch 
ureg <- with(TMB.data ,ifelse(creg==1,1,0))
ur <- colSums(ureg*cpue,na.rm=TRUE)
urp <- colSums(ureg*TMB.data$testf,na.rm=TRUE)

# Restricted Catch 
ureg <- with(TMB.data ,ifelse(creg==2,1,0))
rst <- colSums(ureg*cpue,na.rm=TRUE)
rp <- colSums(ureg*TMB.data $testf,na.rm=TRUE)

# Monofilament Catch 
ureg <- with(TMB.data ,ifelse(creg==3|creg==5,1,0))
nst <- colSums(ureg*cpue,na.rm=TRUE)
np <- colSums(ureg*TMB.data$testf,na.rm=TRUE)


TMB.data$nyear <- with(TMB.data,lyear-fyear+1)
TMB.data$wl_esc <- TMB.data$w_esc[1:2,]
TMB.data$wu_esc <- TMB.data$w_esc[3:8,]
TMB.data$al_esc <- TMB.data$a_esc[1:3,]
TMB.data$au_esc <- TMB.data$a_esc[4:14,]
TMB.data$testp <-  as.matrix(t(cbind(urp,rp,np)))
TMB.data$cpue <-  as.matrix(t(cbind(ur,rst,nst)))
TMB.data$fage <- 4
TMB.data$lage <- 7
TMB.data$efN <- ifelse(TMB.data$efn_h<100,25,100)
names(TMB.data)[names(TMB.data)=='efn_h'] <- 'efN_h'
names(TMB.data)[names(TMB.data)=='efn_e'] <- 'efN_e'
names(TMB.data)[names(TMB.data)=='age_h'] <- 'age_p_h'
names(TMB.data)[names(TMB.data)=='age_e'] <- 'age_p_e'
TMB.data$D <- 5
TMB.data$nbyear <- with(TMB.data, nyear+lage-fage)
TMB.data$nage <- with(TMB.data, lage-fage+1)

makedata(TMB.data,'TMB.dat')


parameters <- list(log_trun=rep(13.5,TMB.data$nyear), log_lwesc=rep(5.0,2), log_uwesc=rep(4.0,6), 
log_laesc = rep(4.0,3), log_uaesc=rep(4.0,11), log_rlwesc=rep(-1.0,2), log_rlaesc = rep(-1.0,3), 
log_ruwesc = rep(-1.0,6), log_ruaesc =rep(-1.0,11), log_q = rep(-10,3), log_rq = rep(-1.0,3),
log_fl = rep(1.0,TMB.data$nyear), log_fu=rep(2.0,TMB.data$nyear), mq=0.8, sq=1.0, 
lnalpha = 2.0, beta = 1.0, phi=0.0, lnR0 = 12.5,ln_R = with(TMB.data, rep(12.0,nyear+lage-fage)),
ln_R_Sd = -1, mk_y = with(TMB.data, rep(2.5,nyear+lage-fage)), me_y = with(TMB.data, rep(5.0,nyear+lage-fage)),
mk_Sd=1, me_Sd=1,N_Sd=1)


#'###########################################################################
# Run the ADMB ----
#'###########################################################################	
# Compile model
setwd(TMBdir)
compile(TMB.model,flags="-O0 -g",DLLFLAGS="",libtmb=FALSE)

dyn.load(dynlib("Chinook_reconst_SS_SR_TMB"))
model <- MakeADFun(TMB.data, parameters, DLL="Chinook_reconst_SS_SR_TMB", silent=FALSE) 
fit <- nlminb(model$par, model$fn, model$gr, control = list(eval.max = 10000, iter.max=10000))


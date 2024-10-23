library(R2admb)

library(coda)

library(tidyr)

library(plyr)

library(dplyr)

setwd("C:/Users/sdlarson/OneDrive - State of Alaska/Desktop/ADMB") #working directory where .tpl is found

setup_admb("C:/Program Files (x86)/ADMB") #Where ADMB .bin files are located

setup_admb("S:/REG3/Kuskokwim Research/Data/Run Reconstruction/Chinook - Kuskokwim River/admb/utilities/mingw64") #where g++ compiler is located



#================================================================
# Option 1 for running ADMB Model if you already have the .dat data file. 
#================================================================


compile_admb('ADFG_Rmodel_toADMB_pooled_ln', safe=FALSE, re=FALSE,verbose=FALSE, 
             
             admb_errors= c("stop","warn","ignore"))  #compiles program

run_admb('ADFG_Rmodel_toADMB_pooled_ln') #runs .tpl

admb <- read_admb('ADFG_Rmodel_toADMB_pooled_ln')  #reads output 

AIC = AIC(admb,k=2)  

coef(admb)

summary(admb)

vcov(admb)

logLik(admb)

deviance(admb)

stdEr(admb)

report$`#N`




#================================================================
# Option 2 for running ADMB Model: This version produces that data 
# file from the King Data .CSV 
#================================================================


# 1. Set up libraries and source functions
#===============================================================
library(R2admb)
library(RColorBrewer)
#source("KuskokwimFunctions.R")
#================================================================
# 2. Read in ADFG Monitoring data
#================================================================
#setwd("C:/Users/njsmith/Desktop/ADMB Model Code2") #working directory where .tpl is found
kusko.data <- read.csv("KingData.csv", h=T)   ####  CSV file that was used for R Model. 
#================================================================
#================================================================
# 3. Arrange data for model input
#================================================================
#================================================================
#########################################################################################
#  3.1 Testfishery: Estimate run proporion of 1976-1983                                #   
#########################################################################################


# Extract testfish data
testf<-kusko.data[substr(names(kusko.data),1,3)=='rpw']
# combine week 8, 9 and 10 and drop 
testf[,8] <- testf[,8]+testf[,9]+testf[,10]
testf <- testf[,-(9:10)]
# Replace NA to mean proporion for each week 
for (i in 1:dim(testf)[2]) {
  testf[is.na(testf[i]),i] <- colMeans(testf,na.rm=T)[i]
}
ny2<-length(kusko.data[,1])
#########################################################################################
#  3.2  Rearange fishing effort and harvest data catch 0 to NA                          #   
#########################################################################################
# Extract weekly commercial effort data 
ceff <-kusko.data[substr(names(kusko.data),1,3)=='cew']
# combine week 8, 9  and drop 
ceff[,6] <- ceff[,6]+ceff[,7]
ceff <- ceff[,-7]
# replace 0 to NA
ceff[ceff == 0] <- NA
# Extract weekly commercial catch data
ccat <-kusko.data[substr(names(kusko.data),1,3)=='chw']
# combine week 8, 9  and drop 
ccat[,6] <- ccat[,6]+ccat[,7]
ccat <- ccat[,-7]
# replace 0 to NA
#ccat[ccat == 0&creg==0] <- NA
# Extract weekly commercial est data
creg <-kusko.data[substr(names(kusko.data),1,3)=='cfw']
# combine week 8, 9  and drop 
creg[,6] <- pmax(creg[,6],creg[,7])
creg <- creg[,-7]
#kusko.data[kusko.data2==0]<-0
#########################################################################################
#  3.3  Recalculate Inriver data                           
#########################################################################################
# Extract Inriver data 
inr <-kusko.data[substr(names(kusko.data),1,3)=='In.']
# Calculate CV
inr$cv <- inr$In.river.sd/inr$In.river

#########################################################################################
#  3.4  Calculate Others                                                
#########################################################################################

tcatch <- rowSums(kusko.data[substr(names(kusko.data),1,2)=='H.'],dims = 1,na.rm=T)
# Extract escapement data
esc <- kusko.data[substr(names(kusko.data),1,2)=='w.'|substr(names(kusko.data),1,2)=='a.']
t.esc <- kusko.data$In.river - tcatch
# Calculate observed minimum escapement
minesc <- rowSums(esc, na.rm=T, dims = 1)
# Calculate observed minimum run
minrun <- rowSums(cbind(tcatch,esc), na.rm=T, dims = 1)

ny <- length(kusko.data[,1])

#########################################################################################
#  3.5 Construct dataset used for likelihood modeling                                   
#########################################################################################

kusko.like.data <- as.matrix(cbind(tcatch,inr,esc,testf[3:8],ccat,ceff,creg))


#########################################################################################
#  4.1 Set Initial value and boundaries                                                 #   
#########################################################################################
# Initial starting point
#init <- c(rep(log(250000),ny),rep(5,6),rep(4,14),rep(-10,3),rep(2,6),rep(2,14))	
# Lower bounds
#lb <-  c(log(minrun),rep(2,6), rep(3,14),rep(-14,3),rep(-3,6),rep(-3,14))
# Upper bounds 
#ub <-  c(rep(log(500000),ny),rep(7,6),rep(8,14),rep(-5,3),rep(5,6),rep(5,14))

#======================================================================================================

#  SWITCH TO ADMB IMPLEMENTATION

#======================================================================================================

creg[creg==4]<-0
creg[creg==5]<-3
#=========================================================
# 6.1 Generate ADMB data file for ADFG data
#=========================================================
adfgdat<-list(
  
  nyear = length(kusko.data[,1]), # number of years of data
  nweek = 6,            # number of harvest weeks?
  nweir = 6,            # number of weirs
  nair = 14,            # number of tributaries with aerial surveys
  
  testf = testf[,3:8],  # proportion of run by week and year
  
  ceff = ceff,          # matrix of commercial effort (weeks in rows, years in columns)
  ccat = ccat,          # matrix of commercial catch (weeks in rows, years in columns)
  creg = creg,          # matrix of commercial harvest regulation indices
  
  inriv = inr$In.river,  # In-river run estimates
  inriv_sd = inr$In.river.sd, # SD of in-river run estimates
  
  tcatch = tcatch,      # Total harvest
  
  esc_w = esc[,1:6],    # Escapement indices from weirs
  esc_a = esc[,7:20],   # Escapement indices from aerial surveys
  
  minesc = minesc,      # minimum escapement (observed at weirs plus aerial surveys)
  minrun = log(minrun), # minimum run size (observed harvest plus minesc)
  ubrun = rep(log(500000),ny)      # upper bound for total run estimates
)

#=====================================================================
# Change NAs back to Zeros in dataset for use in ADMB
#=====================================================================
for(i in 1:15)
{
  adfgdat[[i]][is.na(adfgdat[[i]])]<-0
}


#write_dat("ADFG_Rmodel_toADMB_pooled_ln.dat", adfgdat)


#=====================================================================
# Run ADMB
#=====================================================================

fn<-"ADFG_Rmodel_toADMB_pooled_ln"

start.parms<-list(
  log_trun = rep(12.5,ny2),
  log_wesc = rep(5,6),
  log_aesc = rep(4,14),
  log_q = rep(-11,2),
  log_cvw = 1,
  log_cva = 1,
  log_cvq = 1 )

aout<-do_admb(fn=fn,data=adfgdat,mcmc=F,params=start.parms,
              verbose=T,run.opts = run.control(clean_files = "none",read_files = F))

adparmsout<-c(list(fn = fn),read_pars(fn))

admb <- read_admb('ADFG_Rmodel_toADMB_pooled_ln') 
summary(admb)











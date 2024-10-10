#'==============================================================================
#   Kuskokwim Chinook Salmon Run Reconstruction R-code TMB Implementation 
#
#   Written by  H. Hamaazaki 
#   Date  04/30/2024         
#   # Requirement 
#   # Rtools istalled
#	#Libarary RTMB
#'==============================================================================

#'==============================================================================
#  1.0  Initialize working Environment ----                                        
#'==============================================================================
rm(list=ls(all=TRUE))
library(RTMB)
# set this Year 
this.year <- 2023
# Set Base directry 
Base_dir <- file.path(getwd(),this.year)
# Specify TMB directry 
Model_dir <- file.path(Base_dir,'Model')
# Specify Data directory 
Data_dir <- file.path(Base_dir,'Data')
# Specify Ouptput directory
Out_dir <- file.path(Base_dir,'Outputs')
# Install ADMB read write R source code:
source(file.path(Model_dir,'R_functions','ADMBtoR.r'))
# Install RTMB model
source(file.path(Model_dir,'TMB','Chinook_RR_Model_default.r'))
source(file.path(Model_dir,'R_functions','TMB_functions.r'))
# Set input data file name 
datafn <- file.path(Data_dir,paste0('Chinook_RR_data_',this.year,'.dat'))

create.dataset <- FALSE

if(isTRUE(create.dataset)){
data_file.1 <- 'Kusko_Chinook_RR_Input_2023.csv'
# Filenane for age data  
data_file.2 <- 'Kusko_Chinook_RR_Age_2023.csv'
# Filenane for age data  
data_file.3 <- 'Kusko_Chinook_data_lookup.csv'
source(file.path(Model_dir,'R_functions','Create_data.r'))
source(file.path(Model_dir,'R_functions','Create_data_plot.r'))
source(file.path(Model_dir,'R_functions','Create_age_data.r'))
source(file.path(Model_dir,'R_functions','Create_list_data.r'))
makedata(dat,datafn)
TMB.data <- dat
}else{
#  Read ADMB data  to list 
TMB.data <- datatoR(datafn) 
}
nyear <- with(TMB.data, lyear-fyear+1)
nweir <- dim(TMB.data$w_esc)[1]-1
naerial <- dim(TMB.data$a_esc)[1]
#'------------------------------------------------------------------------------
#  Set Pamateter inital, upper, and lower bounds                                     
#'------------------------------------------------------------------------------
parameters <- list(log_trun=rep(12.5,nyear), log_wesc=rep(5.0,nweir), log_aesc = rep(4.0,naerial), 
log_cvw = 0.5, log_cva=0.5,log_q = rep(-11,2), log_cvq = 0.5)
p.lower <- c(log_trun=TMB.data$minrun, log_wesc=rep(0.0,nweir), log_aesc = rep(0.0,naerial), 
log_cvw = -10.0, log_cva=-10.0,log_q = rep(-12,2), log_cvq = -10.0)
p.upper <- c(log_trun=rep(13.5,nyear), log_wesc=rep(7.0,nweir), log_aesc = rep(7.0,naerial),log_cvw = 1.0, log_cva=1.0, 
log_q = rep(-9,2), log_cvq = 1.0)

#'------------------------------------------------------------------------------
# Run the RTMB  ----
#'------------------------------------------------------------------------------
fit <- runRTMB(TMB.data,Chinook_RR_RTMB_default,parameters,p.upper,p.lower) 
# Compile model
rep <- sdreport(fit)
par1 <-as.data.frame(summary(rep))
names(par1)[2] <- 'Std.Error'
write.csv(par1,file.path(Out_dir,'Chinook_RR_output.csv'))
#'------------------------------------------------------------------------------
# Extract RTMB  Estimates ----
#'------------------------------------------------------------------------------
t.run <- par1[substr(row.names(par1),1,5)=='t_run',]
t.esc <- par1[substr(row.names(par1),1,3)=='esc',]
w.esc <- par1[substr(row.names(par1),1,8)=='log_wesc',]
a.esc <- par1[substr(row.names(par1),1,8)=='log_aesc',]
q.cpue <- par1[substr(row.names(par1),1,5)=='log_q',]
cv <- par1[substr(row.names(par1),1,6)=='log_cv',]

rept <- fit$report()


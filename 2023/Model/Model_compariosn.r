#'==============================================================================
#  Kuskokwim_Chinook_RR_RTMB_Model_Diagnoes.R
#  This program reads RTMB model ouptputs and plot 
#   Written by  H. Hamaazaki 
#   Date  10/25/2024         
#   # Requirement 
#   # Install Rtool
#	# Install Packages TMB, RTMB
#'==============================================================================
#'==============================================================================
#  1.0  Initialize working Environment ----                                        
#'==============================================================================
rm(list=ls(all=TRUE))
library(RTMB)
# Set Base dirctory (May not be needed when running via Rstudio
Base_dir <- file.path('C:','Projects','Kuskokwim_River','Chinook_reconst','Kuskokwim_Chinook')
# set this Year 
this.year <- 2023
Base_dir <- file.path(Base_dir, this.year)
#Base_dir <- file.path('.', this.year)
# Set Base directry 
# Specify TMB directry 
Model_dir <- file.path(Base_dir,'Model')
# Specify Data directory 
Data_dir <- file.path(Base_dir,'Data')
# Specify Ouptput directory
Out_dir <- file.path(Base_dir,'Outputs')
# Install ADMB read write R source code:
source(file.path(Model_dir,'R_functions','ADMBtoR.r'))
source(file.path(Model_dir,'R_functions','TMB_functions.r')) # Needed to run RTMB
#'------------------------------------------------------------------------------
#'  Set Input data (Update every year)
#'------------------------------------------------------------------------------
# RR Input data 
rr_data <- 'Kusko_Chinook_RR_Input_2023.csv'  
# Used to make model 
look_up_data <- 'Kusko_Chinook_lookup.csv'
# Model output 
models <- c('Chinook_RR_output_1.csv', 'Chinook_RR_output_3.csv')
nmodel <- length(models)
# Read to R file 
kusko.data <- read.csv(file.path(Data_dir,rr_data),header=TRUE, na.string='')
lookup <- read.csv(file.path(Data_dir,look_up_data),header=TRUE, na.string='')
output <- list()
for(i in nmodel){
output[[i]] <- read.csv(file.path(Out_dir,models[[i]]),header=TRUE, na.string='')
}


#'------------------------------------------------------------------------------
#'------------------------------------------------------------------------------
# Extract RTMB  Estimates ----
#'------------------------------------------------------------------------------
year <- kusko.data$Year
nyear <- length(year)

t.run <- output[substr(output$X,1,5)=='t_run',]
t.esc <- output[substr(output$X,1,5)=='t_esc',]
log.trun <- output[substr(output$X,1,8)=='log_trun',]
w.esc <- output[substr(output$X,1,8)=='log_wesc',]
a.esc <- output[substr(output$X,1,8)=='log_aesc',]
q.cpue <- output[substr(output$X,1,5)=='log_q',]
cv <- output[substr(output$X,1,6)=='log_cv',]
#rept <- fit$report()
qs <- output[substr(output$X,1,2)=='qs',1]
mq <- output[substr(output$X,1,2)=='mq',1]

fl <- output[substr(output$X,1,6)=='log_fl',1]
fu <- output[substr(output$X,1,6)=='log_fu',1]

t.Run <- t.run$Estimate
Run_uci <- with(log.trun, exp(Estimate+2*Std.Error))
Run_lci <- with(log.trun, exp(Estimate-2*Std.Error))
esc <- t.esc$Estimate
esc_uci <- with(t.esc, exp(log(Estimate)+2*(Std.Error/Estimate)))
esc_lci <- with(t.esc, exp(log(Estimate)-2*(Std.Error/Estimate)))

############################################################################
# 4.1 Plot Run, Escapement, 
############################################################################
u <- 1000
##### Set graphing parameter ########################################################################
windows(record=TRUE)
par(mfrow=c(1,1),mar = c(4, 4, 1, 1))  
plot(year,t.Run/u,type='o',ljoin=1,ylab = 'Total run x 1000',pch=16,yaxt='n', xaxt='n',ylim=c(0,550))
    arrows(year,y0= Run_uci/u,y1=Run_lci/u,code=0)
#: yaxt = 'n'  do not plot the y-axis 
axis(side=2, at = seq(0,550,by= 50),las=2) 
axis(side=1, at = with(dat,seq(fyear,lyear,by= 5))) 
points(year,with(kusko.data,In.river/u),pch=16,col='red')
arrows(year,y0=with(kusko.data,In.river+2*In.river.sd)/u,y1=with(kusko.data,In.river-2*In.river.sd)/u,code=0,col='red',lwd=2)
arrows(year,y0=esc_uci/u,y1=esc_lci/u,code=0, col=1,lty=5)
lines(year,esc/u,type='o',ljoin=1,pch=21,yaxt='n', xaxt='n',col=1,bg='white',lty=5)
legend('topright',legend=c('Run','Escapement','In River Run estimate'), lty=c(1,5,1),pch=c(16,21,16),col=c(1,1,2), bg='white', bty ='n')



############################################################################
# 4.1.1 Plot Upper Escapement & MR
############################################################################
# Plot Run 
plot(year,(rept$up_esc)/1000,type='o',ljoin=1,ylab = 'Upper Escapement x 1000',pch=16, xaxt='n',yaxt='n',ylim=c(0,500), main='Kuskokwim Chinook')
# Plot 95% CI
axis(side=2, at = seq(0,500,by= 50),las=2) 
axis(side=1, at = seq(1975,1975+ny,by= 5)) 
arrows(year,y0=exp(log(kusko.data$mr)+2*kusko.data$mr.sd/kusko.data$mr)/1000,y1=exp(log(kusko.data$mr)-2*kusko.data$mr.sd/kusko.data$mr)/1000,code=0,col='red',lwd=2)
points(year,kusko.data$mr/1000,pch=16,col='red')
legend('topright',legend=c('Upper Escapement','Upper MR'), lty=c(1,1),pch=c(16,16),col=c(1,2), bg='white', bty ='n')

############################################################################
# 4.1.2 Plot Lower Run and Sonar 
############################################################################
# Plot Run 
plot(year,(rept$low_run)/1000,type='o',ljoin=1,ylab = 'Lower Run x 1000',pch=16, xaxt='n',yaxt='n',ylim=c(0,500), main='Kuskokwim Chinook')
# Plot 95% CI
axis(side=2, at = seq(0,500,by= 50),las=2) 
axis(side=1, at = seq(1975,1975+ny,by= 5)) 
# Plot Sonar 
arrows(year,y0=exp(log(kusko.data$Sonar)+2*kusko.data$Sonar.sd/kusko.data$Sonar)/1000,y1=exp(log(kusko.data$Sonar)-2*kusko.data$Sonar.sd/kusko.data$Sonar)/1000,code=0,col='blue',lwd=2)
points(year,kusko.data$Sonar/1000,pch=16,col='blue')
legend('topright',legend=c('Lower Run','Sonar '), lty=c(1,1),pch=c(16,16),col=c(1,2), bg='white', bty ='n')

############################################################################
# 4.2 Plot Harvest rate
############################################################################
plot(year,(P.exp),type='o',ljoin=1,ylab = 'Run Harvest Rate',pch=16,ylim=c(0,1),main='Kuskokwim Chinook Harvest Rate')


############################################################################
# 4.3 Run / Escapement Dignoses 
############################################################################
# Get Weir and Aerial parameters
weir <- c(exp(output[which(output$name=='log_lwesc'),3]),exp(output[which(output$name=='log_uwesc'),3]))
aerial <- c(exp(output[which(output$name=='log_laesc'),3]),exp(output[which(output$name=='log_uaesc'),3]))
sq <- c((output[which(output$name=='sq'),3]))
mq <- c((output[which(output$name=='mq'),3]))
par(mfrow=c(5,5),mar = c(1.75,1.5,1.5,1.75),oma = c(3,3,3,3),cex=0.6) 
# Lower MR
x <-kusko.data$tmr
y <-t.Run
maxx <- max(max(x,na.rm=TRUE)/1000,max(y,na.rm=TRUE)/1000)
plot(year,y/1000,ylab='',type='l', main='Lower MR',ylim=c(0,maxx))
points(year,x/1000,pch=19,col='red')
RMSE <- round(sqrt(sum((log(x)-log(y))^2,na.rm=T)/sum(!is.na(x))),2)
tex <- paste('RMSE ',RMSE,sep = "")
legend('topleft',tex,bty='n',xjust=0,inset=-0.05)

# Upper MR
x <-kusko.data$mr
y <-mq*rept$Escapement
maxx <- max(max(x,na.rm=TRUE)/1000,max(y,na.rm=TRUE)/1000)
plot(year,y/1000,ylab='',type='l', main='Upper MR',ylim=c(0,maxx))
points(year,x/1000,pch=19,col='red')
RMSE <- round(sqrt(sum((log(x)-log(y))^2,na.rm=T)/sum(!is.na(x))),2)
tex <- paste('RMSE ',RMSE,sep = "")
legend('topleft',tex,bty='n',xjust=0,inset=-0.05)

#Sonar
x <-c(kusko.data$Sonar[1:39]/sq[1],kusko.data$Sonar[-c(1:39)]/sq[2])
y <-c(rept$low_run)
maxx <- max(max(x,na.rm=TRUE)/1000,max(y,na.rm=TRUE)/1000)
plot(year,y/1000,ylab='',type='l', main='Sonar',ylim=c(0,maxx))
points(year,x/1000,pch=19,col='red')
RMSE <- round(sqrt(sum((log(x)-log(y))^2,na.rm=T)/sum(!is.na(x))),2)
tex <- paste('RMSE ',RMSE,sep = "")
legend('topleft',tex,bty='n',xjust=0,inset=-0.05)


W.name <- c('Kwethluk Weir','Tuluksak Weir','Aniak Weir','George Weir','Kogrukluk Weir','Tatlawiksuk Weir','Takotna Weir',
 'Pitka Weir')
 
A.name <- c('Kwethluk Aerial','Tuluksak Aerial','Kisaralik Aerial','Salmon (Aniak) Aerial','Kipchuk Aerial','Aniak Aerial',
'Holokuk Aeriai','Oskawalik Aerial','Holitna Aerial','Cheeneetnuk Aerial','Gagaryah Aerial','Pitka Aerial',
'Bear Aerial','Salmon (Pitka) Aerial')
w.esc.n <- names(kusko.data)[substr(names(kusko.data),1,2)=='w.']
a.esc.n <- names(kusko.data)[substr(names(kusko.data),1,2)=='a.']

for (i in 1:8){
x <-t(kusko.data[,w.esc.n[i]])
y <-rept$Escapement/weir[i]

maxx <- max(max(x,na.rm=TRUE),max(y,na.rm=TRUE))
plot(year,y,ylab='',type='l', main=W.name[i],ylim=c(0,maxx))
points(year,x,col='red', pch = 19, cex=1.5)
RMSE <- round(sqrt(sum((log(x)-log(y))^2,na.rm=T)/sum(!is.na(x))),2)
tex <- paste('RMSE ',RMSE,sep = "")
legend('topleft',tex,bty='n',xjust=0,inset=-0.05)
}

for (i in 1:14){
x <-t(kusko.data[a.esc.n[i]])
y <-c(rept$Escapement/aerial[i])
 
 
maxx <- max(max(x,na.rm=TRUE),max(y,na.rm=TRUE))
plot(year,y,ylab='',type='l', main=A.name[i], ylim=c(0,maxx))
points(year,x,col='red', pch = 19, cex=1.5)
RMSE <- round(sqrt(sum((log(x)-log(y))^2,na.rm=T)/sum(!is.na(x))),2)
tex <- paste('RMSE ',RMSE,sep = "")
legend('topleft',tex,bty='n',xjust=0,inset=-0.05)
}

########  Add Texts  #######################################################
mtext("Observed vs. Predicted Counts", side = 3, line = 0, outer = TRUE)
mtext('Predicted', side = 2, line = 1, outer = TRUE)
mtext("Observed", side = 1, line = 1, outer = TRUE)


############################################################################
# 4.4 Fit to Fishery Data 
############################################################################
###### Create Expected Weekly ##############################################
q <- c(exp(output[which(output$name=='log_q'),3]))
  cpue <- as.matrix(ccat/ceff)
# Observed Run proporion adjusted cpue 
  cpue_N <- as.matrix(ccat/ceff/testf[,3:8])
############################################################################
#   Calculate likelihoood for unrestricted 
############################################################################
par(mfrow=c(2,2),mar = c(1.75,1.5,1.5,1.75),oma = c(3,3,3,3),cex=0.6) 
# Unrestricted Mesh 
# Extract all mesh regulation year/week 
unr <- creg
# Keep unrestricted mesh regulation year/week 1: indicate unrestricted period
unr[unr != 1] <- NA
# Keep only cpue of Unrestricted 
unr.N <- cpue*unr
unr.p <- testf[,3:8]*unr
# Rmove all NA 
x <- rowSums(unr.N,na.rm =TRUE)/rowSums(unr.p,na.rm =TRUE)
x[is.nan(x)] <- NA
y <- q[1]*t.Run
maxx <- max(x,y,na.rm=TRUE)
z <-cbind(x,y)
plot(z,xlab = '', ylab ='',xlim=c(0,maxx),ylim=c(0,maxx),main='Unrestricted')
RMSE <- round(sqrt(sum((log(x)-log(y))^2,na.rm=TRUE)/sum(!is.na(x))),2)
tex <- paste('RMSE ',RMSE,sep = "")
legend('topleft',tex,bty='n')
abline(0,1)

# Restricted Mesh 
# Extract all mesh regulation year/week 
unr <- creg
# Keep unrestricted mesh regulation year/week 
unr[unr != 2] <- NA
unr.N <- cpue*unr
unr.p <- testf[,3:8]*unr
# Rmove all NA 
x <- rowSums(unr.N,na.rm =TRUE)/rowSums(unr.p,na.rm =TRUE)
x[is.nan(x)] <- NA
y <- q[2]*t.Run
maxx <- max(x,y,na.rm=TRUE)
z <-cbind(x,y)
plot(z,xlab = '', ylab ='',xlim=c(0,maxx),ylim=c(0,maxx),main='Restricted')
RMSE <- round(sqrt(sum((log(x)-log(y))^2,na.rm=TRUE)/sum(!is.na(x))),2)
tex <- paste('RMSE ',RMSE,sep = "")
legend('topleft',tex,bty='n')
abline(0,1)

# Extract all mesh regulation year/week 
unr <- creg
# Keep unrestricted mesh regulation year/week 1: indicate unrestricted period
unr[(unr != 3)&(unr != 5)] <- NA
unr.N <- cpue*unr
unr.p <- testf[,3:8]*unr
# Rmove all NA 
x <- rowSums(unr.N,na.rm =TRUE)/rowSums(unr.p,na.rm =TRUE)
x[is.nan(x)] <- NA
y <- q[3]*t.Run
maxx <- max(x,y,na.rm=TRUE)
z <-cbind(x,y)
plot(z,xlab = '', ylab ='',xlim=c(0,maxx),ylim=c(0,maxx),main='Restricted Mono')
RMSE <- round(sqrt(sum((log(x)-log(y))^2,na.rm=TRUE)/sum(!is.na(x))),2)
tex <- paste('RMSE ',RMSE,sep = "")
legend('topleft',tex,bty='n')
abline(0,1)

########  Add Texts  #######################################################
mtext("Observed vs. Predicted CPUE", side = 3, line = 0, outer = TRUE)
mtext('Predicted', side = 2, line = 1, outer = TRUE)
mtext("Observed", side = 1, line = 1, outer = TRUE)

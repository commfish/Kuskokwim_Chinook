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
library(RTMB)
#library(TMBdebug)
rm(list=ls(all=TRUE))
# Set Base directry 
Base_dir <- file.path('C:','Projects','Kuskokwim_River','Chinook_reconst','Hamachan')
# Specify TMB directry 
TMBdir <- file.path(Base_dir,'TMB','Model')
# Specify Data directory 
Data_dir <- file.path(Base_dir,'2022')
# Install ADMB read write R source code:
source(file.path(Base_dir,'ADMBtoR.r'))
# Specify TMB Source file name 
TMB.model <- 'Chinook_reconst_TMB.cpp'
# 
data_file <- 'Kusko_Chinook_2022.dat'
datafn <- file.path(Data_dir,data_file)
# Create data set?
create.dataset <- FALSE
if(isTRUE(create.dataset)){
data_file1 <- 'Kusko_RR_Input_Oct_2022.csv'
# Filenane for age data  
data_file2 <- 'Kusko_RR_Age_Oct_2022.csv'
# Filenane for age data  
data_file3 <- 'Kusko_data_lookup.csv'
source(file.path(Base_dir,'Kusko_Create_data.r'))
makedata(dat,datafn)
TMB.data <- dat
}else{
#  Read ADMB data  to list 
TMB.data <- datatoR(datafn) 
}



TMB.data$nyear <- with(TMB.data,lyear-fyear+1)
TMB.data$wl_esc <- TMB.data$w_esc[1:2,]
TMB.data$wu_esc <- TMB.data$w_esc[3:8,]
TMB.data$al_esc <- TMB.data$a_esc[1:3,]
TMB.data$au_esc <- TMB.data$a_esc[4:14,]


parameters <- list(log_trun=rep(13.5,TMB.data$nyear), log_lwesc=rep(5.0,2), log_uwesc=rep(4.0,6), 
log_laesc = rep(4.0,3), log_uaesc=rep(4.0,11), log_rlwesc=rep(-1.0,2), log_rlaesc = rep(-1.0,3), 
log_ruwesc = rep(-1.0,6), log_ruaesc =rep(-1.0,11), log_q = rep(-10,3), log_rq = rep(-1.0,3),
log_fl = rep(1.0,TMB.data$nyear), log_fu=rep(2.0,TMB.data$nyear), mq=0.2, sq=1.0)

p.lower <- c(log_trun=rep(10,TMB.data$nyear), log_lwesc=rep(1.0,2), log_uwesc=rep(1.0,6), 
log_laesc = rep(1.0,3), log_uaesc=rep(1.0,11), log_rlwesc=rep(-10.0,2), log_rlaesc = rep(-10.0,3), 
log_ruwesc = rep(-10.0,6), log_ruaesc =rep(-10.0,11), log_q = rep(-12,3), log_rq = rep(-10.0,3),
log_fl = rep(0.0,TMB.data$nyear), log_fu=rep(0.0,TMB.data$nyear), mq=0.0, sq=0.5)

p.upper <- c(log_trun=rep(14.5,TMB.data$nyear), log_lwesc=rep(7.0,2), log_uwesc=rep(7.0,6), 
log_laesc = rep(7.0,3), log_uaesc=rep(7.0,11), log_rlwesc=rep(1.0,2), log_rlaesc = rep(1.0,3), 
log_ruwesc = rep(1.0,6), log_ruaesc =rep(1.0,11), log_q = rep(-9.0,3), log_rq = rep(2.0,3),
log_fl = rep(10.0,TMB.data$nyear), log_fu=rep(10.0,TMB.data$nyear), mq=1.0, sq=2.0)

map <- list(log_rlwesc=rep(factor(NA),2), log_rlaesc = rep(factor(NA),3), 
log_ruwesc = rep(factor(NA),6), log_ruaesc =rep(factor(NA),11), log_rq = rep(factor(NA),3))

phase2Init <- c( best[names(best) %in% c('log_trun','log_lwesc','log_uwesc','log_laesc','log_uaesc')],log_rlwesc=rep(1.0,2), log_rlaesc = rep(1.0,3), 
log_ruwesc = rep(1.0,6), log_ruaesc =rep(1.0,11),best[names(best) %in% c('log_q')],log_rq = rep(2.0,3), best[names(best) %in% c('log_fl','log_fu','mq','sq')])

length(phase2Init)



#'###########################################################################
# Run the TMB ----
#'###########################################################################	
# Compile model
setwd(TMBdir)
compile(TMB.model)

dyn.load(dynlib("Chinook_reconst_TMB"))
model <- MakeADFun(TMB.data, parameters, DLL="Chinook_reconst_TMB", silent=FALSE) 
fit <- nlminb(model$par, model$fn, model$gr, control = list(eval.max = 10000, iter.max=10000))


model <- MakeADFun(TMB.data, parameters, DLL="Chinook_reconst_TMB", map=map,silent=FALSE) 
fit <- nlminb(model$par, model$fn, model$gr, control = list(eval.max = 10000, iter.max=10000))
best <- model$env$last.par.best
phase2Init <- c( best[names(best) %in% c('log_trun','log_lwesc','log_uwesc','log_laesc','log_uaesc')],log_rlwesc=rep(1.0,2), log_rlaesc = rep(1.0,3), 
log_ruwesc = rep(1.0,6), log_ruaesc =rep(1.0,11),best[names(best) %in% c('log_q')],log_rq = rep(2.0,3), best[names(best) %in% c('log_fl','log_fu','mq','sq')])
model <- MakeADFun(TMB.data, parameters, DLL="Chinook_reconst_TMB", silent=FALSE) 
fit <- nlminb(phase2Init, model$fn, model$gr, control = list(eval.max = 10000, iter.max=10000),
lower = p.lower, upper = p.upper)
best <- model$env$last.par.best


dyn.load(dynlib("Chinook_reconst_TMB"))
model <- MakeADFun(TMB.data, parameters, DLL="Chinook_reconst_TMB", silent=FALSE) 
fit <- nlminb(model$par, model$fn, model$gr, control = list(eval.max = 10000, iter.max=10000),lower = p.lower, upper = p.upper)
rep <- sdreport(model)
print(summary(rep,p.value=T))
par1 <-summary(rep)
write.csv(par1,'temp.csv')



print(best)
print(rep)

names(fit)

# How many iterations did it take to converge?
fit$iterations

# What was my objective function value at our MLE
fit$objective

# Do a likelihood profile
prof <- tmbprofile(model,"log_trun",trace=F)
plot(prof)





#'###########################################################################
# Extract Data ---
#'###########################################################################
dat <- rapply(TMB.data, f=function(x) ifelse(x==0,NA,x), how="replace" )

year <- with(TMB.data,c(fyear:lyear))

# Extract N
t.Run <- exp(par1[which(rownames(par1)=='log_trun'),1])
# Extract S
Esc <- par1[which(rownames(par1)=='t_esc'),]
P.exp <- (t.Run-Esc$value)/t.Run

Run_uci <- exp(par1[which(rownames(par1)=='log_trun'),1]+2*par1[which(rownames(par1)=='log_trun'),2])
Run_lci <- exp(par1[which(rownames(par1)=='log_trun'),1]-2*par1[which(rownames(par1)=='log_trun'),2]) 
Esc_uci <- exp(log(Esc[,1])+2*(Esc[,2]/Esc[,1]))
Esc_lci <- exp(log(Esc[,1])-2*(Esc[,2]/Esc[,1]))
# Create Total run, escapement, lci, uci
run.table <- data.frame(cbind(t.Run,Run_lci,Run_uci,Esc$value,Esc_lci,Esc_uci))
# Write table to CSV file
#write.csv(run.table,file=paste(admodeldir,'run_table.csv',sep=''),na='')

dat$mr_tr
dat$mrsd_tr
tmr <- with(dat,data.frame(mr_tr,mrsd_tr))


###########################################################################
#    Standard Model Diagnostics 
############################################################################
##### Set graphing parameter ###############################################
windows(record=TRUE)
par(mfrow=c(1,1),mar = c(4, 4, 1, 1))  

############################################################################
# 4.1 Plot Run, Escapement, 
############################################################################
#windows(record=TRUE)
#par(mfrow=c(1,1),mar = c(4, 4, 1, 1)) 
ny <- length(year) 
# Plot Run 
plot(year,(t.Run)/1000,type='o',ljoin=1,ylab = 'Drainage Run x 1000',pch=16, xaxt='n',yaxt='n',ylim=c(0,500), main='Kuskokwim Chinook Upper Escapement')
# Plot 95% CI
arrows(year,y0=Run_uci/1000,y1=Run_lci/1000,code=0)
axis(side=2, at = seq(0,500,by= 50),las=2) 
axis(side=1, at = seq(1975,1975+ny,by= 5)) 
arrows(year,y0=exp(with(dat,log(mr_tr)+2*mrsd_tr/mr_tr))/1000,y1=exp(with(dat,log(mr_tr)-2*mrsd_tr/mr_tr))/1000,code=0,col='red',lwd=2)
points(year,dat$mr_tr/1000,pch=16,col='red')
# Plot Sonar 
sq <- par1[which(rownames(par1)=='sq'),1]
sh <- dat$sonar/sq+dat$h_low
arrows(year,y0=exp(with(dat,log(sh)+2*sonar_sd/sonar))/1000,y1=exp(with(dat,log(sh)-2*sonar_sd/sonar))/1000,code=0,col='blue',lwd=2)
points(year,sh/1000,pch=16,col='blue')

# Plot Escapement
arrows(year,y0=Esc_uci/1000,y1=Esc_lci/1000,code=0, col= 4)
lines(year,(Esc[,1])/1000,type='o',ljoin=1,pch=16, col=4)
legend('topright',legend=c('Run','Escapement','Lower MR'), lty=c(1,1,1),pch=c(16,16,16),col=c(1,4,2), bg='white', bty ='n')
mq <- par1[which(rownames(par1)=='mq'),1]
points(year,dat$mr/(1-mq)/1000,pch=16,col='red')

############################################################################
# 4.1.1 Plot Upper Escapement & MR
############################################################################
# Plot Run 
mq <- par1[which(rownames(par1)=='mq'),1]
plot(year,(mq*Esc[,1])/1000,type='o',ljoin=1,ylab = 'Upper Escapement x 1000',pch=16, xaxt='n',yaxt='n',ylim=c(0,500), main='Kuskokwim Chinook')
# Plot 95% CI
axis(side=2, at = seq(0,500,by= 50),las=2) 
axis(side=1, at = seq(1975,1975+ny,by= 5)) 
arrows(year,y0=exp(log(dat$mr)+2*dat$mr.sd/dat$mr)/1000,y1=exp(log(dat$mr)-2*dat$mr.sd/dat$mr)/1000,code=0,col='red',lwd=2)
points(year,dat$mr/1000,pch=16,col='red')
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
arrows(year,y0=exp(log(dat$Sonar)+2*dat$Sonar.sd/dat$Sonar)/1000,y1=exp(log(dat$Sonar)-2*dat$Sonar.sd/dat$Sonar)/1000,code=0,col='blue',lwd=2)
points(year,dat$Sonar/1000,pch=16,col='blue')
legend('topright',legend=c('Lower Run','Sonar '), lty=c(1,1),pch=c(16,16),col=c(1,2), bg='white', bty ='n')

############################################################################
# 4.2 Plot Harvest rate
############################################################################
plot(year,(P.exp),type='o',ljoin=1,ylab = 'Run Harvest Rate',pch=16,ylim=c(0,1),main='Kuskokwim Chinook Harvest Rate')


############################################################################
# 4.3 Run / Escapement Dignoses 
############################################################################
# Get Weir and Aerial parameters
weir <- c(exp(par1[which(rownames(par1)=='log_lwesc'),1]),exp(par1[which(rownames(par1)=='log_uwesc'),1]))
aerial <- c(exp(par1[which(rownames(par1)=='log_laesc'),1]),exp(par1[which(rownames(par1)=='log_uaesc'),1]))

par(mfrow=c(5,5),mar = c(1.75,1.5,1.5,1.75),oma = c(3,3,3,3),cex=0.6) 
# Lower MR
x <-dat$mr_tr
y <-t.Run
maxx <- max(max(x,na.rm=TRUE)/1000,max(y,na.rm=TRUE)/1000)
plot(year,y/1000,ylab='',type='l', main='Lower MR',ylim=c(0,maxx))
points(year,x/1000,pch=19,col='red')
RMSE <- round(sqrt(sum((log(x)-log(y))^2,na.rm=T)/sum(!is.na(x))),2)
tex <- paste('RMSE ',RMSE,sep = "")
legend('topleft',tex,bty='n',xjust=0,inset=-0.05)

# Upper MR
x <-dat$mr
y <-(1-mq)*Esc[,1]
maxx <- max(max(x,na.rm=TRUE)/1000,max(y,na.rm=TRUE)/1000)
plot(year,y/1000,ylab='',type='l', main='Upper MR',ylim=c(0,maxx))
points(year,x/1000,pch=19,col='red')
RMSE <- round(sqrt(sum((log(x)-log(y))^2,na.rm=T)/sum(!is.na(x))),2)
tex <- paste('RMSE ',RMSE,sep = "")
legend('topleft',tex,bty='n',xjust=0,inset=-0.05)

#Sonar
x <-c(dat$sonar[1:39]/sq,dat$sonar[-c(1:39)]/sq)
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
w.esc.n <- dat$w_esc
a.esc.n <- dat$a_esc


for (i in 1:8){
x <-w.esc.n[i,]
if(i <3){
y <-mq*Esc[,1]/(weir[i])
} else {
y <-(1-mq)*Esc[,1]/(weir[i])
}
maxx <- max(max(x,na.rm=TRUE),max(y,na.rm=TRUE))
plot(year,y,ylab='',type='l', main=W.name[i],ylim=c(0,maxx))
points(year,x,col='red', pch = 19, cex=1.5)
RMSE <- round(sqrt(sum((log(x)-log(y))^2,na.rm=T)/sum(!is.na(x))),2)
tex <- paste('RMSE ',RMSE,sep = "")
legend('topleft',tex,bty='n',xjust=0,inset=-0.05)
}

for (i in 1:14){
x <-a.esc.n[i,]
if(i <4){
y <-mq*Esc[,1]/(aerial[i])
} else {
y <-(1-mq)*Esc[,1]/(aerial[i])
}
 
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
q <- c(exp(par1[which(par1$name=='log_q'),3]))
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

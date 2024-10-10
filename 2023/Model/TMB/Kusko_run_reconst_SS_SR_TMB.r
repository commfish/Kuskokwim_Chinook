#'==============================================================================
#   Kuskokwim Chinook Salmon State-Space Model Run Reconstruction R-code 
#   TMB Implementation 
#
#   Written by  H. Hamaazaki 
#   Date  07/03/2023         
#     
#	This 
#
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
# Set TMB directry 
TMBdir <- file.path(Base_dir,'Model')
# Set Data directory 
Data_dir <- file.path(Base_dir,'Data')
# Install ADMB read write R source code:
source(file.path('C:','Projects','Statistics_Worksheets','R_functions','ADMBtoR.r'))


# Set TMB Source file name 
TMB.model <- 'Chinook_reconst_SS_SR_TMB.cpp'


# Specify the name of accompanying program file 
datafn <- file.path(Data_dir,'Kusko_Chinook_22_new_v.dat')
#  Read ADMB data  to list 
TMB.data <- datatoR(datafn) 



TMB.data$nyear <- with(TMB.data,lyear-fyear+1)
TMB.data$wl_esc <- TMB.data$w_esc[1:2,]
TMB.data$wu_esc <- TMB.data$w_esc[3:8,]
TMB.data$al_esc <- TMB.data$a_esc[1:3,]
TMB.data$au_esc <- TMB.data$a_esc[4:14,]
TMB.data$efN <- ifelse(TMB.data$efn_h<100,25,100)
names(TMB.data)[names(TMB.data)=='efn_h'] <- 'efN_h'
names(TMB.data)[names(TMB.data)=='efn_e'] <- 'efN_e'
names(TMB.data)[names(TMB.data)=='age_h'] <- 'age_p_h'
names(TMB.data)[names(TMB.data)=='age_e'] <- 'age_p_e'
TMB.data$D <- 5
TMB.data$nbyear <- with(TMB.data, nyear+lage-fage)
TMB.data$nage <- with(TMB.data, lage-fage+1)




parameters <- list(log_trun=rep(13.5,TMB.data$nyear), log_lwesc=rep(5.0,2), log_uwesc=rep(4.0,6), 
log_laesc = rep(4.0,3), log_uaesc=rep(4.0,11), log_rlwesc=rep(-1.0,2), log_rlaesc = rep(-1.0,3), 
log_ruwesc = rep(-1.0,6), log_ruaesc =rep(-1.0,11), log_q = rep(-10,3), log_rq = rep(-1.0,3),
log_fl = rep(1.0,TMB.data$nyear), log_fu=rep(2.0,TMB.data$nyear), mq=0.8, sq=1.0, 
lnalpha = 2.0, beta = 1.0, phi=0.0, lnR0 = 12.5,ln_R = with(TMB.data, rep(12.0,lage)), omega = with(TMB.data, rep(0,nyear-fage)),ln_R_Sd = -1, mk_y = with(TMB.data, rep(2.5,nyear+lage-fage)), me_y = with(TMB.data, rep(5.0,nyear+lage-fage)))

random <- c('mk_y','me_y','omega')

map <- list(log_rlwesc=rep(factor(NA),2), log_rlaesc = rep(factor(NA),3), 
log_ruwesc = rep(factor(NA),6), log_ruaesc =rep(factor(NA),11), log_rq = rep(factor(NA),3))


p.lower <- c(log_trun=rep(10,TMB.data$nyear), log_lwesc=rep(1.0,2), log_uwesc=rep(1.0,6), 
log_laesc = rep(1.0,3), log_uaesc=rep(1.0,11), log_rlwesc=rep(-10.0,2), log_rlaesc = rep(-10.0,3), 
log_ruwesc = rep(-10.0,6), log_ruaesc =rep(-10.0,11), log_q = rep(-12,3), log_rq = rep(-10.0,3),
log_fl = rep(0.0,TMB.data$nyear), log_fu=rep(0.0,,TMB.data$nyear), mq=0.5, sq=0.5)

p.upper <- c(log_trun=rep(14.5,TMB.data$nyear), log_lwesc=rep(7.0,2), log_uwesc=rep(7.0,6), 
log_laesc = rep(7.0,3), log_uaesc=rep(7.0,11), log_rlwesc=rep(1.0,2), log_rlaesc = rep(1.0,3), 
log_ruwesc = rep(1.0,6), log_ruaesc =rep(1.0,11), log_q = rep(-9.0,3), log_rq = rep(2.0,3),
log_fl = rep(10.0,TMB.data$nyear), log_fu=rep(10.0,TMB.data$nyear), mq=1.0, sq=2.0)

map <- list(log_rlwesc=rep(factor(NA),2), log_rlaesc = rep(factor(NA),3), 
log_ruwesc = rep(factor(NA),6), log_ruaesc =rep(factor(NA),11), log_rq = rep(factor(NA),3))


phase2Init <- c( best[names(best) %in% c('log_trun','log_lwesc','log_uwesc','log_laesc','log_uaesc')],log_rlwesc=rep(1.0,2), log_rlaesc = rep(1.0,3), 
log_ruwesc = rep(1.0,6), log_ruaesc =rep(1.0,11),best[names(best) %in% c('log_q')],log_rq = rep(2.0,3), best[names(best) %in% c('log_fl','log_fu','mq','sq')])

length(phase2Init)

fixwinpath <- function() {
  PATH <- Sys.getenv("PATH")
  PATH <- paste0(R.home(), "/bin/x64;", PATH)
  PATH <- paste0("c:/rtools42/mingw64/bin;", PATH)
  Sys.setenv(PATH=PATH)
}

#'###########################################################################
# Run the ADMB ----
#'###########################################################################	
# Compile model
parameters <- list(log_trun=rep(13.5,TMB.data$nyear), log_lwesc=rep(5.0,2), log_uwesc=rep(4.0,6), 
log_laesc = rep(4.0,3), log_uaesc=rep(4.0,11), log_rlwesc=rep(-1.0,2), log_rlaesc = rep(-1.0,3), 
log_ruwesc = rep(-1.0,6), log_ruaesc =rep(-1.0,11), log_q = rep(-10,3), log_rq = rep(-1.0,3),
log_fl = rep(1.0,TMB.data$nyear), log_fu=rep(2.0,TMB.data$nyear), mq=0.8, sq=1.0, 
lnalpha = 2.0, beta = 1.0, phi=0.0,ln_R = with(TMB.data, rep(12.0,lage)), omega = with(TMB.data, rep(0,nyear-fage)), mk_y = with(TMB.data, rep(2.5,nyear+lage-fage)), me_y = with(TMB.data, rep(5.0,nyear+lage-fage)))

setwd(TMBdir)
compile(TMB.model,flags="-O0 -g",DLLFLAGS="",random=random,libtmb=FALSE)
dyn.load(dynlib("Chinook_reconst_SS_SR_TMB"))
model <- MakeADFun(TMB.data, parameters, DLL="Chinook_reconst_SS_SR_TMB", silent=FALSE) 
fit <- nlminb(model$par, model$fn, model$gr, control = list(eval.max = 10000, iter.max=10000))

rep <- sdreport(model)
print(summary(rep))
par1 <-summary(rep)
rept <- model$report()


t.Run <- exp(par1[which(rownames(par1)=='log_trun'),1])
# Extract S
Esc <- par1[which(rownames(par1)=='S'),]
N <- par1[which(rownames(par1)=='N'),]
year <- with(TMB.data,c(fyear:lyear))
plot(year,(t.Run)/1000,type='o',ljoin=1,ylab = 'Drainage Run x 1000',pch=16, xaxt='n',ylim=c(0,500), main='Kuskokwim Chinook')
lines(year,N[,1]/1000,lty=2,col=3)

best <- model$env$last.par.best
phase2Init <- c( best[names(best) %in% c('log_trun','log_lwesc','log_uwesc','log_laesc','log_uaesc')],log_rlwesc=rep(1.0,2), log_rlaesc = rep(1.0,3), 
log_ruwesc = rep(1.0,6), log_ruaesc =rep(1.0,11),best[names(best) %in% c('log_q')],log_rq = rep(2.0,3), best[names(best) %in% c('log_fl','log_fu','mq','sq')])
model <- MakeADFun(TMB.data, parameters, DLL="Chinook_reconst_TMB", silent=FALSE) 
fit <- nlminb(phase2Init, model$fn, model$gr, control = list(eval.max = 10000, iter.max=10000),
lower = p.lower, upper = p.upper)
best <- model$env$last.par.best



model <- MakeADFun(TMB.data, parameters, DLL="Chinook_reconst_TMB", silent=FALSE) 
fit <- nlminb(model$par, model$fn, model$gr, control = list(eval.max = 10000, iter.max=10000))


print(rep)

names(fit)

# How many iterations did it take to converge?
fit$iterations

# What was my objective function value at our MLE
fit$objective

# Do a likelihood profile
prof <- tmbprofile(model,"log_trun",trace=F)
plot(prof)




print(summary(rep,p.value=T))
par1 <-summary(rep)
write.csv(par1,'temp.csv)



#'###########################################################################
# Extract Data ---
#'###########################################################################

year <- with(TMB.data,c(fyear:lyear))
plot(year,model$report()$t_run)
lines(year, (model$report()$N))
N_pa <- model$report()$N_pa
N_pa_ob <- model$report()$N_pa_ob

plot(year,N_pa_ob[,2])
lines(year,N_pa[,2],pch=19)


# Extract N
t.Run <- exp(par1[which(rownames(par1)=='log_trun'),1])
# Extract S
Esc <- par1[which(rownames(par1)=='S'),]
N <- par1[which(rownames(par1)=='N'),]

P.exp <- (t.Run-Esc[1])/t.Run

Run_uci <- exp(par1[which(rownames(par1)=='log_trun'),1]+2*par1[which(rownames(par1)=='log_trun'),2])
Run_lci <- exp(par1[which(rownames(par1)=='log_trun'),1]-2*par1[which(rownames(par1)=='log_trun'),2]) 
Esc_uci <- exp(log(Esc[,1])+2*(Esc[,2]/Esc[,1]))
Esc_lci <- exp(log(Esc[,1])-2*(Esc[,2]/Esc[,1]))
# Create Total run, escapement, lci, uci
run.table <- data.frame(cbind(t.Run,Run_lci,Run_uci,Esc$value,Esc_lci,Esc_uci))
# Write table to CSV file
#write.csv(run.table,file=paste(admodeldir,'run_table.csv',sep=''),na='')
mq <- par1[which(rownames(par1)=='mq'),1]

dat <- rapply(TMB.data, f=function(x) ifelse(x==0,NA,x), how="replace" )
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
plot(year,(t.Run)/1000,type='o',ljoin=1,ylab = 'Drainage Run x 1000',pch=16, xaxt='n',yaxt='n',ylim=c(0,500), main='Kuskokwim Chinook')
lines(year,N[,1]/1000,lty=2,col=3)
# Plot 95% CI
arrows(year,y0=Run_uci/1000,y1=Run_lci/1000,code=0)
axis(side=2, at = seq(0,500,by= 50),las=2) 
axis(side=1, at = seq(1975,1975+ny,by= 5)) 
arrows(year,y0=exp(log(dat$mr_tr)+2*dat$mrsd_tr/dat$mr_tr)/1000,y1=exp(log(dat$mr_tr)-2*dat$mrsd_tr/dat$mr_tr)/1000,code=0,col='red',lwd=2)
points(year,dat$mr_tr/1000,pch=16,col='red')
# Plot Sonar 
arrows(year,y0=exp(log(dat$sonar)+2*dat$sonar_sd/dat$sonar)/1000,y1=exp(log(dat$sonar)-2*dat$sonar_sd/dat$sonar)/1000,code=0,col='blue',lwd=2)
points(year,dat$sonar/1000,pch=16,col='blue')
# Plot Upper MR
arrows(year,y0=exp(log(dat$mr)+2*dat$mrsd/dat$mr)/mq/1000,y1=exp(log(dat$mr)-2*dat$mrsd/dat$mr)/mq/1000,code=0,col='red',lwd=2)
points(year,dat$mr/mq/1000,pch=16,col='red')
# Plot Escapement
arrows(year,y0=Esc_uci/1000,y1=Esc_lci/1000,code=0, col= 4)
lines(year,(Esc[,1])/1000,type='o',ljoin=1,pch=16, col=4)
legend('topright',legend=c('Run','Escapement','Lower MR','Sonar', ), lty=c(1,1,1),pch=c(16,16,16),col=c(1,4,2), bg='white', bty ='n')


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
legend('topright',legend=c('Escapement','Upper MR'), lty=c(1,1),pch=c(16,16),col=c(1,2), bg='white', bty ='n')

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
sq <- c((par1[which(rownames(par1)=='sq'),1]))
mq <- c((par1[which(rownames(par1)=='mq'),1]))
par(mfrow=c(5,5),mar = c(1.75,1.5,1.5,1.75),oma = c(3,3,3,3),cex=0.6) 
# Lower MR
x <-dat$tmr
y <-t.Run
maxx <- max(max(x,na.rm=TRUE)/1000,max(y,na.rm=TRUE)/1000)
plot(year,y/1000,ylab='',type='l', main='Lower MR',ylim=c(0,maxx))
points(year,x/1000,pch=19,col='red')
RMSE <- round(sqrt(sum((log(x)-log(y))^2,na.rm=T)/sum(!is.na(x))),2)
tex <- paste('RMSE ',RMSE,sep = "")
legend('topleft',tex,bty='n',xjust=0,inset=-0.05)

# Upper MR
x <-dat$mr
y <-mq*Esc[,1]
maxx <- max(max(x,na.rm=TRUE)/1000,max(y,na.rm=TRUE)/1000)
plot(year,y/1000,ylab='',type='l', main='Upper MR',ylim=c(0,maxx))
points(year,x/1000,pch=19,col='red')
RMSE <- round(sqrt(sum((log(x)-log(y))^2,na.rm=T)/sum(!is.na(x))),2)
tex <- paste('RMSE ',RMSE,sep = "")
legend('topleft',tex,bty='n',xjust=0,inset=-0.05)

#Sonar
x <-c(dat$Sonar[1:39]/sq[1],dat$Sonar[-c(1:39)]/sq[2])
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
w.esc.n <- names(dat)[substr(names(dat),1,2)=='w.']
a.esc.n <- names(dat)[substr(names(dat),1,2)=='a.']

for (i in 1:8){
x <-t(dat[,w.esc.n[i]])
y <-Esc[,1]/weir[i]

maxx <- max(max(x,na.rm=TRUE),max(y,na.rm=TRUE))
plot(year,y,ylab='',type='l', main=W.name[i],ylim=c(0,maxx))
points(year,x,col='red', pch = 19, cex=1.5)
RMSE <- round(sqrt(sum((log(x)-log(y))^2,na.rm=T)/sum(!is.na(x))),2)
tex <- paste('RMSE ',RMSE,sep = "")
legend('topleft',tex,bty='n',xjust=0,inset=-0.05)
}

for (i in 1:14){
x <-t(dat[a.esc.n[i]])
y <-c(Esc[,1]/aerial[i])
 
 
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


############################################################################
#  3.2  Calculate SR Parameters 
############################################################################
# Ricker alpha 
	alpha <- exp(lnalpha$value)
# Ricker alpha.c 	
	lnalpha.c <- lnalpha$value+var(Rec.dev$value)/(1-phi$value^2)/2
# Ricker beta
	beta <- sbeta$value/100000
# Ricker Smsy	
	Smsy <- lnalpha$value*(0.5-0.07*lnalpha$value)/beta
# Ricker Seq
	Seq <- lnalpha$value/beta	

############################################################################
#  4.5 Plot Run age composition
############################################################################
par(mfrow=c(1,1),mar = c(4, 4, 1, 1),oma=c(1,1,1,1))  
# Plot Run 
plot(year,(rept$N_ta[1,])/1000,type='o',ljoin=1,ylab = 'Drainage Run x 1000',pch=16, xaxt='n',yaxt='n',ylim=c(0,500), main='Kuskokwim Chinook')
lines(year,(rept$N_ta[2,])/1000,type='o',ljoin=1,pch=16, col=2)
lines(year,(rept$N_ta[3,])/1000,type='o',ljoin=1,pch=16, col=3)
lines(year,(rept$N_ta[4,])/1000,type='o',ljoin=1,pch=16, col=4)
axis(side=2, at = seq(0,500,by= 50),las=2) 
axis(side=1, at = seq(1975,1975+ny,by= 5)) 
# Plot Legend
legend('topright',legend=c('Age 4','Age 5','Age 6','Age 7'), lty=1,pch=16,col=c(1,2,3,4), bg='white', bty ='n')

par(mfrow=c(1,1),mar = c(4, 4, 1, 1),oma=c(1,1,1,1))  
# Plot Run 
plot(year,(rept$N_pa[,1]),type='o',ljoin=1,ylab = 'Drainage Run x 1000',pch=16, xaxt='n',yaxt='n',ylim=c(0,1), main='Kuskokwim Chinook')
lines(year,(rept$N_pa[,2]),type='o',ljoin=1,pch=16, col=2)
lines(year,(rept$N_pa[,3]),type='o',ljoin=1,pch=16, col=3)
lines(year,(rept$N_pa[,4]),type='o',ljoin=1,pch=16, col=4)
axis(side=2, at = seq(0,1,by= .2),las=2) 
axis(side=1, at = seq(1975,1975+ny,by= 5)) 
# Plot Legend
legend('topright',legend=c('Age 4','Age 5','Age 6','Age 7'), lty=1,pch=16,col=c(1,2,3,4), bg='white', bty ='n')


# Extract Harvest age proptition
p.est.age<- rept$N_pa
# Extract Escapment  age proptition
p.obs.age<- rept$N_pa_ob

par(mfrow=c(7,7),mar = c(1, 1, 1, 1),oma=c(4,4,2,1),mgp=c(3,.3,0))
for (i in 1:ny){
x <- plot(c(4:7),p.obs.age[i,],ylim=c(0,0.8),tck = -0.05,pch=19, xaxt='n',bty='l')
axis(side=1,tck = -0.05, cex=0.8, at = seq(4,7)) 
lines(c(4:7),p.est.age[i,],lty=2)
text(6.5,0.8,paste(year[i]),cex=1.0)
}
mtext("Age Composition: observed vs predicted", side=3, outer=TRUE)
mtext('Age proporion', side = 2, line = 1, outer = TRUE)
mtext("Age", side = 1, line = 1, outer = TRUE)




############################################################################
#  Figure for bubble plot
############################################################################
# Create commercial length vs predicted length freq figure 
bubbleplot<-function(pdata,obdata){
bdata <- as.matrix(obdata - pdata)
fbdata <- as.vector(bdata)
classes <- sort(rep(4:7,length(year)))
years <- (rep(year,4))
test2 <- data.frame(cbind(years,classes,fbdata))
test2$radius <- 2*sqrt(abs(test2$fbdata)/pi) 
symbols(test2$years,test2$classes,circles=test2$radius, bg=ifelse(test2$fbdata<0,'white','black'), xlim=c(min(years),max(years)),tck = -0.02, inches = FALSE,xaxt='n',yaxt='n',xlab = "year", ylab="age")
axis(side=1, at = seq(fyear,lyear,by= 5)) 
# Plot x bar
axis(side=2, at = c(4,5,6,7)) 
}
par(mfrow=c(2,1),mar = c(4, 4, 4, 4))  
bubbleplot(t(p.est.age),t(p.obs.age))



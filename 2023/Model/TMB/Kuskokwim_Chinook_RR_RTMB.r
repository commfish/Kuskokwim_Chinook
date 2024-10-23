#'==============================================================================
#  Kuskokwim_Chinook_RR_RTMB.R 
#  Kuskokwim Chinook Salmon Run Reconstruction R-code RTMB Implementation 
#  This model runs the RTMB version of the Kuskokwim Run reconstruction model 
#   Written by  H. Hamaazaki 
#   Date  04/30/2024         
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
Base_dir <- file.path(Base_dir,this.year)
# Set Base directry 
# Specify TMB directry 
Model_dir <- file.path(Base_dir,'Model')
# Specify Data directory 
Data_dir <- file.path(Base_dir,'Data')
# Specify Ouptput directory
Out_dir <- file.path(Base_dir,'Outputs')
# Install ADMB read write R source code:
source(file.path(Model_dir,'R_functions','ADMBtoR.r'))

#'------------------------------------------------------------------------------
#'  Input data  
#'------------------------------------------------------------------------------

rr_data <- 'Kusko_Chinook_RR_Input_2023.csv'
age_data <- 'Kusko_Chinook_RR_age_2023.csv'
look_up_data <- 'Kusko_Chinook_lookup.csv'
kusko.data <- read.csv(file.path(Data_dir,rr_data),header=TRUE, na.string='')
age_data <- read.csv(file.path(Data_dir,age_data),header=TRUE, na.string='')
look_up_data <- read.csv(file.path(Data_dir,look_up_data),header=TRUE, na.string='')
#'------------------------------------------------------------------------------

create.dataset <- TRUE
if(isTRUE(create.dataset)){
source(file.path(Model_dir,'R_functions','Create_RR_data.r'))
dat <- make_RR_Model_data(kusko.data)
# Filename for age data  
#source(file.path(Model_dir,'R_functions','Create_data_plot.r'))
#source(file.path(Model_dir,'R_functions','Create_age_data.r'))
makedata(dat,datafn)
TMB.data <- dat
}else{
#  Read ADMB data  to list 
TMB.data <- datatoR(datafn) 
}
nyear <- with(TMB.data, lyear-fyear+1)
nweir <- dim(TMB.data$w_esc)[1]
naerial <- dim(TMB.data$a_esc)[1]


#'Run default model ------------------------------------------------------------
#'------------------------------------------------------------------------------
#  Set Pamateter inital, upper, and lower bounds                                     
#'------------------------------------------------------------------------------
# Orignal model + Salmon(Pitka) Weir 
parameters <- list(log_trun=rep(12.5,nyear), log_wesc=rep(5.0,nweir), 
log_aesc = rep(4.0,naerial), log_cvw =  -0.7, log_cva= -0.7,log_q = rep(-11,2), 
log_cvq = -0.7)
p.lower <- c(log_trun=log(TMB.data$minrun), log_wesc=rep(0.0,nweir), 
log_aesc = rep(0.0,naerial), log_cvw = -10.0, log_cva=-10.0,log_q = rep(-12,2), 
log_cvq = -10.0)
p.upper <- c(log_trun=rep(13.5,nyear), log_wesc=rep(7.0,nweir), 
log_aesc = rep(7.0,naerial),log_cvw = 1.0, log_cva=1.0, log_q = rep(-9,2),log_cvq = 1.0)

#'------------------------------------------------------------------------------
# Run the RTMB  ----
#'------------------------------------------------------------------------------
source(file.path(Model_dir,'TMB','Chinook_RR_Model_RTMB_default.r'))
fit <- runRTMB(TMB.data,Chinook_RR_Model_RTMB_default,parameters,p.upper,p.lower) 
#'End default model Run --------------------------------------------------------
# Compile model
rep <- sdreport(fit)
par1 <-as.data.frame(summary(rep))
names(par1)[2] <- 'Std.Error'
write.csv(par1,file.path(Out_dir,'Chinook_RR_output_1.csv'))


#'Run default_rev model --------------------------------------------------------
#'------------------------------------------------------------------------------
#  Set Pamateter inital, upper, and lower bounds                                     
#'------------------------------------------------------------------------------
parameters <- list(log_trun=rep(12.5,nyear), log_wesc=rep(5.0,nweir), 
log_aesc = rep(4.0,naerial), log_cvw = -0.7, log_cva=-0.7,log_q = rep(-11,2), 
log_cvq = -0.7,log_fl = rep(1.0, nyear), log_fu = rep(2.0,nyear))

p.lower <- c(log_trun=log(TMB.data$minrun), log_wesc=rep(0.0,nweir), 
log_aesc = rep(0.0,naerial), log_cvw = -10.0, log_cva=-10.0,log_q = rep(-12,2), 
log_cvq = -10.0,log_fl = rep(-10.0, nyear), log_fu = rep(-10.0,nyear))

p.upper <- c(log_trun=rep(13.5,nyear), log_wesc=rep(7.0,nweir), 
log_aesc = rep(7.0,naerial),log_cvw = 1.0, log_cva=1.0, log_q = rep(-9,2), 
log_cvq = 1.0,log_fl = rep(5.0, nyear), log_fu = rep(5.0,nyear))
#'------------------------------------------------------------------------------
# Run the RTMB  ----
#'------------------------------------------------------------------------------
source(file.path(Model_dir,'TMB','Chinook_RR_Model_RTMB_default_rev.r'))
fit <- runRTMB(TMB.data,Chinook_RR_Model_RTMB_default_rev,parameters,p.upper,p.lower) 
#'End default model Run --------------------------------------------------------
# Compile model
rep <- sdreport(fit)
par1 <-as.data.frame(summary(rep))
names(par1)[2] <- 'Std.Error'
write.csv(par1,file.path(Out_dir,'Chinook_RR_output_2.csv'))

#'Run New model ------------------------------------------------------------
#'------------------------------------------------------------------------------
#  Set Pamateter inital, upper, and lower bounds                                     
#'------------------------------------------------------------------------------
parameters <- list(log_trun=rep(12.5,nyear), log_wesc=rep(5.0,nweir), 
log_aesc = rep(4.0,naerial), log_cvw = 0.5, log_cva=0.5,log_q = rep(-11,2), 
log_cvq = 0.5,mq = 0.25, log_fl = rep(1.0, nyear), log_fu = rep(2.0,nyear))
p.lower <- c(log_trun=TMB.data$minrun, log_wesc=rep(0.0,nweir), 
log_aesc = rep(0.0,naerial), log_cvw = -10.0, log_cva=-10.0,log_q = rep(-12,2), 
log_cvq = -10.0,mq=0.1,log_fl = rep(-10.0, nyear), log_fu = rep(-10.0,nyear))
p.upper <- c(log_trun=rep(13.5,nyear), log_wesc=rep(7.0,nweir), 
log_aesc = rep(7.0,naerial),log_cvw = 1.0, log_cva=1.0, log_q = rep(-9,2), 
log_cvq = 1.0,mq=0.5,log_fl = rep(5.0, nyear), log_fu = rep(5.0,nyear))

#'------------------------------------------------------------------------------
# Run the RTMB  ----
#'------------------------------------------------------------------------------
source(file.path(Model_dir,'TMB','Chinook_RR_Model_RTMB_New.r'))
fit <- runRTMB(TMB.data,Chinook_RR_Model_RTMB_New,parameters,p.upper,p.lower) 
#'End default model Run --------------------------------------------------------

rep <- sdreport(fit)
par1 <-as.data.frame(summary(rep))
names(par1)[2] <- 'Std.Error'
write.csv(par1,file.path(Out_dir,'Chinook_RR_output_3.csv'))



#'------------------------------------------------------------------------------
# Extract RTMB  Estimates ----
#'------------------------------------------------------------------------------
year <- with(dat,c(fyear:lyear))
t.run <- par1[substr(row.names(par1),1,5)=='t_run',]
t.esc <- par1[substr(row.names(par1),1,5)=='t_esc',]
log.trun <- par1[substr(row.names(par1),1,8)=='log_trun',]
w.esc <- par1[substr(row.names(par1),1,8)=='log_wesc',]
a.esc <- par1[substr(row.names(par1),1,8)=='log_aesc',]
q.cpue <- par1[substr(row.names(par1),1,5)=='log_q',]
cv <- par1[substr(row.names(par1),1,6)=='log_cv',]
rept <- fit$report()
qs <- par1[substr(row.names(par1),1,2)=='qs',1]
mq <- par1[substr(row.names(par1),1,2)=='mq',1]

fl <- par1[substr(row.names(par1),1,6)=='log_fl',1]
fu <- par1[substr(row.names(par1),1,6)=='log_fu',1]

t.Run <- t.run$Estimate[(nyear+1):(2*nyear)]
Run_uci <- with(log.trun, exp(Estimate+2*Std.Error))[(nyear+1):(2*nyear)]
Run_lci <- with(log.trun, exp(Estimate-2*Std.Error))[(nyear+1):(2*nyear)]
esc <- t.esc$Estimate[(nyear+1):(2*nyear)]
esc_uci <- with(t.esc, exp(log(Estimate)+2*(Std.Error/Estimate)))[(nyear+1):(2*nyear)]
esc_lci <- with(t.esc, exp(log(Estimate)-2*(Std.Error/Estimate)))[(nyear+1):(2*nyear)]


#####################################################################################################
#    Plot scatter graphs
#####################################################################################################
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


########### Scatter Graph : Weir/Aerial vs. predicted escapement #####################################
par(mfrow=c(3,3),mar = c(1.75,1.5,1.5,1.75),oma = c(3,3,3,3),cex=0.6) 
t.Run <- t.run$Estimate[(1):(nyear)]
x <-t(kusko.data[,'tmr'])
y <-t.Run
maxx <- max(max(x,na.rm=TRUE),max(y,na.rm=TRUE))
plot(year,y,ylab='',type='l', main='Lower MR',ylim=c(0,maxx))
points(year,x,col='red', pch = 19, cex=1.5)
RMSE <- round(sqrt(sum((log(x)-log(y))^2,na.rm=T)/sum(!is.na(x))),2)
tex <- paste('RMSE ',RMSE,sep = "")
legend('topleft',tex,bty='n',xjust=0,inset=-0.05)

x <-t(kusko.data[,'sonar'])
y <-t.Run*(1-exp(-fl))
maxx <- max(max(x,na.rm=TRUE),max(y,na.rm=TRUE))
plot(year,y,ylab='',type='l', main='Bethel Sonar',ylim=c(0,maxx))
points(year,x,col='red', pch = 19, cex=1.5)
RMSE <- round(sqrt(sum((log(x)-log(y))^2,na.rm=T)/sum(!is.na(x))),2)
tex <- paste('RMSE ',RMSE,sep = "")
legend('topleft',tex,bty='n',xjust=0,inset=-0.05)

x <- TMB.data$h_low
y <-t.Run*(exp(-fl))
maxx <- max(max(x,na.rm=TRUE),max(y,na.rm=TRUE))
plot(year,y,ylab='',type='l', main='Lower Catch',ylim=c(0,maxx))
points(year,x,col='red', pch = 19, cex=1.5)
RMSE <- round(sqrt(sum((log(x)-log(y))^2,na.rm=T)/sum(!is.na(x))),2)
tex <- paste('RMSE ',RMSE,sep = "")
legend('topleft',tex,bty='n',xjust=0,inset=-0.05)


x <- TMB.data$h_up
y <-t.Run*(1-exp(-fl))*(exp(-fu))
maxx <- max(max(x,na.rm=TRUE),max(y,na.rm=TRUE))
plot(year,y,ylab='',type='l', main='Upper Catch',ylim=c(0,maxx))
points(year,x,col='red', pch = 19, cex=1.5)
RMSE <- round(sqrt(sum((log(x)-log(y))^2,na.rm=T)/sum(!is.na(x))),2)
tex <- paste('RMSE ',RMSE,sep = "")
legend('topleft',tex,bty='n',xjust=0,inset=-0.05)

x <- t(kusko.data[,'mr'])
y <-(1-mq)*t.Run*(1-exp(-fl))*(1-exp(-fu))
maxx <- max(max(x,na.rm=TRUE),max(y,na.rm=TRUE))
plot(year,y,ylab='',type='l', main='Upper MR',ylim=c(0,maxx))
points(year,x,col='red', pch = 19, cex=1.5)
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
weir <- exp(w.esc$Estimate)
for (i in 1:8){
x <-t(kusko.data[,w.esc.n[i]])
y <-esc/weir[i]
maxx <- max(max(x,na.rm=TRUE),max(y,na.rm=TRUE))
plot(year,y,ylab='',type='l', main=W.name[i],ylim=c(0,maxx))
points(year,x,col='red', pch = 19, cex=1.5)
RMSE <- round(sqrt(sum((log(x)-log(y))^2,na.rm=T)/sum(!is.na(x))),2)
tex <- paste('RMSE ',RMSE,sep = "")
legend('topleft',tex,bty='n',xjust=0,inset=-0.05)
}

aerial <- exp(a.esc$Estimate) 
for (i in 1:14){
x <-t(kusko.data[a.esc.n[i]])
y <-c(esc/aerial[i])

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




###### Create Expected Weekly ##################################################
q1 <- output[which(output$name == 'log_q1'),2]
q2 <- output[which(output$name == 'log_q2'),2]
q3 <- output[which(output$name == 'log_q3'),2]
nweek<-dim(creg)[2]
week.eff <- matrix(0,ncol=nweek,nrow=ny)
for (i in 1:ny){
	for (j in 1:nweek) 
	{
		if(creg[i,j]==1) 
		{
			week.eff[i,j] <- (-log(1-ccat[i,j]/(testf[i,j+2]*output$mean[i]))/q1)
		}         
		if(creg[i,j]==2) 
		{
            week.eff[i,j] <- (-log(1-ccat[i,j]/(testf[i,j+2]*output$mean[i]))/q2)
		}
		if(creg[i,j]==3|creg[i,j]==5) 
		{
            week.eff[i,j] <- (-log(1-ccat[i,j]/(testf[i,j+2]*output$mean[i]))/q3)
		}
	}
}  
windows(record=TRUE)
eff <- as.matrix(ceff)
par(mfrow=c(6,7),mar = c(4, 2, 1, 1))
for (i in 1:ny){
	plot(week.eff[i,],type='l',main=year[i],xlab = "")
	for (j in 1:nweek)
    {
		if(creg[i,j]==1) 
        {
			points(j,eff[i,j],pch=16,col = 'Red')
		}         
		if(creg[i,j]==2) 
		{
			points(j,eff[i,j],pch=16,col = 'blue')
		}
		if(creg[i,j]==3) 
		{
			points(j,eff[i,j],pch=16,col = 'green')
		}
		if(creg[i,j]==5) 
		{
			points(j,eff[i,j],pch=16,col = 'black')
		}
	}
}
################################################################################
################################################################################
################################################################################

############################################################################
#    Standard Model Diagnostics 
############################################################################
##### Set graphing parameter ###############################################
windows(record=TRUE)
par(mfrow=c(1,1),mar = c(4, 4, 1, 1))  

############################################################################
# 4.1 Plot Run, Escapement, 
############################################################################
windows(record=TRUE)
par(mfrow=c(1,1),mar = c(4, 4, 1, 1)) 
ny <- length(year) 
# Plot Run 
plot(year,(t.Run)/1000,type='o',ljoin=1,ylab = 'Drainage Run x 1000',pch=16, xaxt='n',yaxt='n',ylim=c(0,500), main='Kuskokwim Chinook')
# Plot 95% CI
arrows(year,y0=Run_uci/1000,y1=Run_lci/1000,code=0)
axis(side=2, at = seq(0,500,by= 50),las=2) 
axis(side=1, at = seq(1975,1975+ny,by= 5)) 
arrows(year,y0=exp(log(kusko.data$tmr)+2*kusko.data$tmr.sd/kusko.data$tmr)/1000,y1=exp(log(kusko.data$tmr)-2*kusko.data$tmr.sd/kusko.data$tmr)/1000,code=0,col='red',lwd=2)
points(year,kusko.data$tmr/1000,pch=16,col='red')
# Plot Escapement
arrows(year,y0=Esc_uci/1000,y1=Esc_lci/1000,code=0, col= 4)
lines(year,(Esc$value)/1000,type='o',ljoin=1,pch=16, col=4)
# Observed 
lines(year,rept$up_esc/1000,type='o',ljoin=1,pch=21,yaxt='n', xaxt='n',ylim=c(0,500),col=3,bg='white',lty=5)
points(year,kusko.data$tmr/1000,pch=16,col='red')
arrows(year,y0=exp(log(kusko.data$tmr)+2*kusko.data$tmr.sd/kusko.data$tmr)/1000,y1=exp(log(kusko.data$tmr)-2*kusko.data$tmr.sd/kusko.data$tmr)/1000,code=0,col='red',lwd=2)
# mr above Aniak
points(year,kusko.data$mr/1000,pch=16,col=2)
arrows(year,y0=exp(log(kusko.data$mr)+2*kusko.data$mr.sd/kusko.data$mr)/1000,y1=exp(log(kusko.data$mr)-2*kusko.data$mr.sd/kusko.data$mr)/1000,code=0,col=2,lwd=1,lty=2)
# Sonar above Aniak
points(year,(kusko.data$Sonar+kusko.data$H.Sub.l)/1000,pch=16,col=6)
arrows(year,y0=(exp(log(kusko.data$Sonar)+2*kusko.data$Sonar.sd/kusko.data$Sonar)+kusko.data$H.Sub.l)/1000,y1=(exp(log(kusko.data$Sonar)-2*kusko.data$Sonar.sd/kusko.data$Sonar)+kusko.data$H.Sub.l)/1000,code=0,col=6,lwd=1,lty=2)

# Plot Legend

legend('topright',legend=c('Run','Escapement','observed estimate'), lty=c(1,1,1),pch=c(16,16,16),col=c(1,4,2), bg='white', bty ='n')

############################################################################
# 4.2 Plot Harvest rate
############################################################################
plot(year,(P.exp),type='o',ljoin=1,ylab = 'Run Harvest Rate',pch=16,ylim=c(0,1),main='Kuskokwim Chinook Harvest Rate')


############################################################################
# 4.3 Run / Escapement Dignoses 
############################################################################
# Get Weir and Aerial parameters
weir <- c(exp(par1[which(par1$name=='log_lwesc'),3]),exp(par1[which(par1$name=='log_uwesc'),3]))
aerial <- c(exp(par1[which(par1$name=='log_laesc'),3]),exp(par1[which(par1$name=='log_uaesc'),3]))
sq <- c((par1[which(par1$name=='sq'),3]))

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
y <-rept$up_esc
maxx <- max(max(x,na.rm=TRUE)/1000,max(y,na.rm=TRUE)/1000)
plot(year,y/1000,ylab='',type='l', main='Upper MR',ylim=c(0,maxx))
points(year,x/1000,pch=19,col='red')
RMSE <- round(sqrt(sum((log(x)-log(y))^2,na.rm=T)/sum(!is.na(x))),2)
tex <- paste('RMSE ',RMSE,sep = "")
legend('topleft',tex,bty='n',xjust=0,inset=-0.05)

#Sonar
x <-kusko.data$Sonar
y <-c(rept$low_run[1:39]*sq[1],rept$low_run[-c(1:39)]*sq[2])
maxx <- max(max(x,na.rm=TRUE)/1000,max(y,na.rm=TRUE)/1000)
plot(year,y/1000,ylab='',type='l', main='Sonar',ylim=c(0,maxx))
points(year,x/1000,pch=19,col='red')
RMSE <- round(sqrt(sum((log(x)-log(y))^2,na.rm=T)/sum(!is.na(x))),2)
tex <- paste('RMSE ',RMSE,sep = "")
legend('topleft',tex,bty='n',xjust=0,inset=-0.05)


esc.name <- c('Kwethluk Weir','Tuluksak Weir','George Weir','Kogrukluk Weir','Tatlawiksuk Weir','Takotna Weir',
'Kwethluk Aerial','Tuluksak Aerial','Kisaralik Aerial','Salmon (Aniak) Aerial','Kipchuk Aerial','Aniak Aerial',
'Holokuk Aeriai','Oskawalik Aerial','Holitna Aerial','Cheeneetnuk Aerial','Gagaryah Aerial','Pitka Aerial',
'Bear Aerial','Salmon (Pitka) Aerial')
w.esc.n <- names(kusko.data)[substr(names(kusko.data),1,2)=='w.']
a.esc.n <- names(kusko.data)[substr(names(kusko.data),1,2)=='a.']

for (i in 1:6){
x <-t(kusko.data[w.esc.n[i]])
if(i<=2) 
  {
y <-rept$low_esc/weir[i]
  }
else {
y <-rept$up_esc/weir[i]
 } 
maxx <- max(max(x,na.rm=TRUE),max(y,na.rm=TRUE))
plot(year,y,ylab='',type='l', main=esc.name[i],ylim=c(0,maxx))
points(year,x,col='red', pch = 19, cex=1.5)
RMSE <- round(sqrt(sum((log(x)-log(y))^2,na.rm=T)/sum(!is.na(x))),2)
tex <- paste('RMSE ',RMSE,sep = "")
legend('topleft',tex,bty='n',xjust=0,inset=-0.05)
}

for (i in 1:14){
x <-t(kusko.data[a.esc.n[i]])
if(i<=2) 
  {
y <-rept$low_esc/weir[i]/aerial[i]
  } else if(i==3){
y <-rept$low_esc/aerial[i]
 } else if(i < 13){
y <-rept$up_esc/aerial[i]
 } else 
y <-c(rept$up_esc[1:39]/aerial[i],rept$up_esc[-c(1:39)]/aerial[i+2])
 
 
maxx <- max(max(x,na.rm=TRUE),max(y,na.rm=TRUE))
plot(year,y,ylab='',type='l', main=esc.name[i+6], ylim=c(0,maxx))
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

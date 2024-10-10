#'==============================================================================
#Figure - Total Run, Total Escapement, Scalars
#'==============================================================================
#'------------------------------------------------------------------------------
#' Data read 
#'
#'------------------------------------------------------------------------------
rm(list=ls(all=TRUE))
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
source(file.path(Model_dir,'TMB','Model','Chinook_RR_RTMB_default.r'))
source(file.path(Model_dir,'R_functions','TMB_functions.r'))


datafn <- file.path(Data_dir,paste0('Chinook_RR_data_',this.year,'.dat'))
TMB.data <- datatoR(datafn) 
nyear <- with(TMB.data, lyear-fyear+1)
nweir <- dim(TMB.data$w_esc)[1]-1
naerial <- dim(TMB.data$a_esc)[1]

kusko.data <- rapply(TMB.data, f=function(x) ifelse(x==0,NA,x), how="replace" )
par1<-read.csv(file.path(Out_dir,"Chinook_RR_output.csv"), header=TRUE)
t.run <- par1[substr(par1$X,1,5)=='t_run',]
t.esc <- par1[substr(par1$X,1,3)=='esc',]
w.esc <- par1[substr(par1$X,1,8)=='log_wesc',]
a.esc <- par1[substr(par1$X,1,8)=='log_aesc',]
q.cpue <- par1[substr(par1$X,1,5)=='log_q',]
cv <- par1[substr(par1$X,1,6)=='log_cv',]
                        


par(family="serif") #set font to Times New Roman
windows()
year <-with(TMB.data,fyear:lyear)
ny <- length(year) # of run size estimates
ci <- function(data){
    lci <- with(data,exp(log(Estimate)-2*Std.Error/Estimate))
    uci <- with(data,exp(log(Estimate)+2*Std.Error/Estimate))
   out <- data.frame(lci=lci,uci=uci)
   return(out)
}


###############################################################################
#Figure - Total Run, Total Escapement, Scalars
##############################################################################
jpeg("Run and Esc.jpeg", width = 7, height = 4, units = "in", res = 150, pointsize = 10)

t.run.ci <- ci(t.run)
esc.ci <- ci(t.esc)
par(mfrow=c(1,1),mar = c(4, 4, 1, 1)) 
plot(year,t.run$Estimate/1000,type='o',ljoin=1,ylab = 'Total Run x 1,000',pch=16,yaxt='n', xlab="Year",
      xaxt='n',ylim=c(0,550))
arrows(year,y0=t.run.ci$lci/1000,y1=t.run.ci$uci/1000,code=0)
axis(side=2, at = seq(0,550,by= 50),las=2) 
axis(side=1, at = seq(1976,1976+ny,by= 2),las=2) 
points(year,(kusko.data$inriv/1000),pch=16,col="azure4")
arrows(year,y0=(kusko.data$inriv+2*kusko.data$inriv_sd)/1000,y1=(kusko.data$inriv-2*kusko.data$inriv_sd)/1000,code=0,col="azure4",lwd=2)
arrows(year,y0=esc.ci$lci/1000,y1=esc.ci$uci/1000,code=0, col=1,lty=5)
lines(year,t.esc$Estimate/1000,type='o',ljoin=1,pch=21,yaxt='n', xaxt='n',ylim=c(0,550),
      col=1,bg='white',lty=5)
legend('topright',legend=c('Run','Escapement','In River Run Estimate'), 
       lty=c(1,5,1),pch=c(16,21,16),col=c(1,1,"azure4"), bg=
         'white', bty ='n')

dev.off()
###############################################################################
#Figure - CV, number of projects 
##############################################################################
jpeg("CV vs. # of Projects.jpeg", width = 7, height = 4, units = "in", res = 150, pointsize = 10)

par(mar = c(5, 5, 3, 5))
projects<-rowSums(!is.na(kusko.data[,8:27]))#number of projects operated by year. 


lab<-year

bp<-barplot(model.results$CV[1:ny], yaxt="n", xlab="Year", ylab="Coefficient of Variation (CV)",
            axis.lty=1, ylim=c(0,0.30), names.arg=lab, las=2,xaxs = "i",xlim = c(-1, 56.5), cex.axis = 1.4, cex.names=.9)
box()
mtext("Number of Assessment Projects", side = 4, line = 3)

lines(x=bp, y= projects/100, lty=3, lwd=2)
lines(x=bp[-ny], y= rep(mean(model.results$CV[1:ny-1]), ny-1), lty=1, lwd=2) #average cv of dataset

axis(2, at=seq(0,.25,by=.05), labels=paste(100*seq(0,.25,by=.05), "%") )

axis(4, at=seq(0,.25,by=.05), labels=paste(100*seq(0,.25,by=.05)))

dev.off()
###############################################################################
#Figure Scatter Graph : Weir/Aerial vs. predicted escapement #####################################
###############################################################################
jpeg("OBS vs PREDIC by year.jpeg", width = 7, height = 7, units = "in", res = 150)

par(mfrow=c(5,5),mar = c(1.75,1.5,1.5,1.75),oma = c(3,3,3,3)) 
x <-kusko.data$inriv  #  Observed total run  
y <-t.run$Estimate  #  Model predicted 
maxx <- max(max(x,na.rm=TRUE),max(y,na.rm=TRUE))
plot(year,y,ylab='',type='l', main='In River',ylim=c(0,maxx), cex.main =1)
polygon(c(year,rev(year)),c(t.run.ci$uc,rev(t.run.ci$lc)),col=tcol("blue",70),border=NA)
points(year,x,pch=19,col='black')

esc.name <- c('Kwethluk Weir','Tuluksak Weir','George Weir','Kogrukluk Weir','Tatlawiksuk Weir','Takotna Weir',
              'Kwethluk Aerial','Kisaralik Aerial','Tuluksak Aerial','Salmon (Aniak) Aerial','Kipchuk Aerial','Aniak Aerial',
              'Holokuk Aerial','Oskawalik Aerial','Holitna Aerial','Cheeneetnuk Aerial','Gagaryah Aerial','Pitka Aerial',
              'Bear Aerial','Salmon (Pitka) Aerial')
ob.esc <-rbind(kusko.data$w_esc[-8,],kusko.data$a_esc)
esc.pars <- rbind(w.esc,a.esc)
for (i in 1:20){
  x <-ob.esc[i,]
  y <-t.esc$Estimate/as.numeric(exp(esc.pars$Estimate[i]))
  uc <-esc.ci$uc/as.numeric(exp(esc.pars$Estimate[i]))
  lc <-esc.ci$lc/as.numeric(exp(esc.pars$Estimate[i]))
  maxx <- max(max(x,na.rm=TRUE),max(y,na.rm=TRUE))
  plot(year,y,ylab='',type='l', main=esc.name[i],ylim=c(0,maxx), cex.main=1)
  polygon(c(year,rev(year)),c(uc,rev(lc)),col=tcol("blue",70),border=NA)
  points(year,x,pch=19,col='black')
}

mtext("Observed vs. Estimated Counts", side = 3, line = 0, outer = TRUE)
mtext('Abundance', side = 2, line = 1, outer = TRUE)
mtext("Year", side = 1, line = 1, outer = TRUE)


dev.off()

###############################################################################
#Figure - Scatter Graph : Weir/Aerial vs. predicted escapement
###############################################################################
jpeg("OBS vs PREDIC.jpeg", width = 7, height = 7, units = "in", res = 150)

par(mfrow=c(5,5),mar = c(1.75,1.5,1.5,1.75),oma = c(3,3,3,3),cex=0.6) 
x <-kusko.data[,6]
y <-(model.results$Run[1:ny])
maxx <- max(x,na.rm=TRUE)
z <-cbind(x,y)
plot(z,ylab='',col='black', main='In River',xlim=c(0,maxx),ylim=c(0,maxx),cex.main=1)
abline(0,1)
esc.name <- c('Kwethluk Weir','Tuluksak Weir','George Weir','Kogrukluk Weir',
              'Tatlawiksuk Weir','Takotna Weir','Kwethluk Aerial','Kisaralik Aerial',
              'Tuluksak Aerial','Salmon (Aniak) Aerial','Kipchuk Aerial','Aniak Aerial',
              'Holokuk Aerial','Oskawalik Aerial','Holitna Aerial','Cheeneetnuk Aerial',
              'Gagaryah Aerial','Pitka Aerial',
              'Bear Aerial','Salmon (Pitka) Aerial')
for (i in 1:20){
  x <-kusko.data[7+i]
  y <-model.results$Escapement[1:ny]/as.numeric(exp(model.results$par[ny+i]))
  maxx <- max(x,na.rm=TRUE)
  z <-cbind(x,y)
  plot(z,ylab='',col='black', main=esc.name[i], xlim=c(0,maxx),ylim=c(0,maxx),cex.main=1)
  points(z[ny,],col='black', pch = 19, cex=1.5)
  RMSE <- round(sqrt(sum((log(x)-log(y))^2,na.rm=T)/sum(!is.na(x))),2)
  tex <- paste('RMSE ',RMSE,sep = "")
  legend('topleft',tex,bty='n',xjust=0,inset=-0.05)
  abline(0,1)
}
#  Add Texts to Figure
mtext("Observed vs. Estimated Counts", side = 3, line = 0, outer = TRUE)
mtext('Predicted', side = 2, line = 1, outer = TRUE)
mtext("Observed", side = 1, line = 1, outer = TRUE)

dev.off()
###############################################################################
#Figure - Esc. by project overlayed with drainage escapement
###############################################################################

jpeg(" Indv Project total esc.jpeg", width = 7, height = 8, units = "in", res = 150)
install.packages("terra")
library(terra)
esc.name <- c('Kwethluk Weir','Tuluksak Weir','George Weir','Kogrukluk Weir','Tatlawiksuk Weir','Takotna Weir',
              'Kwethluk Aerial','Kisaralik Aerial','Tuluksak Aerial','Salmon (Aniak) Aerial','Kipchuk Aerial','Aniak Aerial',
              'Holokuk Aerial','Oskawalik Aerial','Holitna Aerial','Cheeneetnuk Aerial','Gagaryah Aerial','Pitka Aerial',
              'Bear Aerial','Salmon (Pitka) Aerial')

par(mfrow=c(1,1),mar = c(4, 4, 1, 1))
for (i in 1:20){
  x <-(ob.esc[i,]*as.numeric(exp(esc.pars$Estimate[i])))/1000
  z <- cbind(year,x)
  if(i==1) plot(z[c((ny-2):ny),],pch=19,ylim=c(0,400), xlim=c(2020.5, 2023.5),xaxt="n", xlab="Year", col="azure4",
                ylab='Model Predicted Escapement (x 1,000)')
  points(z[c((ny-2):ny),], pch = 19, cex=1,col="azure4")
  
  #lines(z[c((ny-2):ny),])
  text(z[c((ny-2):ny),], labels=esc.name[i], cex= 0.5, pos=4)
}
points(year[c((ny-2):ny)],t.esc$Estimate[c((ny-2):ny)]/1000,col = 'black', pch = 19, 
       cex=1.5)  
arrows(year[c((ny-2):ny)],y0=esc.ci$uci[c((ny-2):ny)]/1000,
       y1=esc.ci$lci[c((ny-2):ny)]/1000,code=0,col='black',lty=1,lwd=5)
axis(1,at=c(2021,2022,2023))

dev.off()


#install.packages("plotrix")
#library(plotrix)
#jpeg(" Indv Project total esc.jpeg", width = 7, height = 8, units = "in", res = 150)
#library(lattice)
#par(mfrow=c(1,1),mar = c(4, 4, 1, 1))
#for (i in 1:20){
 # x <-(kusko.data[7+i]*as.numeric(exp(model.results$par[ny+i])))/1000
  #z <- cbind(kusko.data$Year,x)
  #if(i==1) xz<-xyplot(z[,2]~z[,1],pch=19,ylim=c(0,350), xlim=c(2015.5, 2018.5),xaxt="n", xlab="Year", col="azure4",
   #               ylab='Model Predicted Escapement (x 1000)', scales=list(x=list( at= c (2016,2017,2018))))
  #jitter(points(z[c((ny-2):ny),], pch = 19, cex=1))
  #lines(z[c((ny-2):ny),])
  #jitter(text(z[c((ny-2):ny),], labels=esc.name[i], cex= 0.5, pos=4))
#}
#points(kusko.data$Year[c((ny-2):ny)],model.results$Escapement[c((ny-2):ny)]/1000,col = 'black', pch = 19, 
 #      cex=1.5)  
#arrows(kusko.data$Year[c((ny-2):ny)],y0=model.results$UCIEsc[c((ny-2):ny)]/1000,
 #      y1=model.results$LCIEsc[c((ny-2):ny)]/1000,code=0,col='black',lty=1,lwd=5)
#axis(1,at=c(2015,2016,2017,2018))

#dev.off()






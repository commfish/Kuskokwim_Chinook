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

kusko.data <- read.csv(file.path(Data_dir,'Kusko_Chinook_RR_Input_2023.csv'),header=TRUE, na.string='')
model.results<-read.csv(file.path(Out_dir,'Kusko_Chinook_RR_output_2023.csv'),header=TRUE, na.string='') #run reconstrucion model results
lookup <- read.csv(file.path(Data_dir,'Kusko_Chinook_lookup.csv'),header=TRUE, na.string='')
###############################################################################
# Data creation 
##############################################################################
Year <- kusko.data$Year
t.run <- model.results[substr(model.results$X,1,5)=='t_run',]
t.esc <- model.results[substr(model.results$X,1,3)=='esc',]
log.trun <- model.results[substr(model.results$X,1,8)=='log_trun',]
w.esc <- model.results[substr(model.results$X,1,8)=='log_wesc',]
a.esc <- model.results[substr(model.results$X,1,8)=='log_aesc',]
q.cpue <- model.results[substr(model.results$X,1,5)=='log_q',]
cv <- model.results[substr(model.results$X,1,6)=='log_cv',]

Run <- t.run$Estimate
UCIRun <- with(log.trun, exp(Estimate+2*Std.Error))
LCIRun <- with(log.trun, exp(Estimate-2*Std.Error))
Esc <- t.esc$Estimate
UCIEsc <- with(t.esc, exp(log(Estimate)+2*(Std.Error/Estimate)))
LCIEsc <- with(t.esc, exp(log(Estimate)-2*(Std.Error/Estimate)))
CV <- t.run$Std.Error/t.run$Estimate

prj.names <- names(kusko.data)[substr(names(kusko.data),1,2) %in% c('a.','w.')]
esc.name <- lookup$Label[lookup$CSV %in% prj.names]
scl <- c(w.esc$Estimate,a.esc$Estimate) 


par(family="serif") #set font to Times New Roman
ny <- length(kusko.data[,1]) # of run size estimates
u <- 1000

###############################################################################
#Figure - Total Run, Total Escapement, Scalars
##############################################################################
jpeg("Run and Esc.jpeg", width = 7, height = 4, units = "in", res = 150, pointsize = 10)

par(mfrow=c(1,1),mar = c(4, 4, 1, 1)) 

plot(Year,Run/u,type='o',ljoin=1,ylab = 'Total Run x 1,000',pch=16,yaxt='n', xlab="Year",
      xaxt='n',ylim=c(0,550))
arrows(Year,y0= UCIRun[1:ny]/u,y1=LCIRun[1:ny]/u,code=0)
axis(side=2, at = seq(0,550,by= 50),las=2) 
axis(side=1, at = seq(min(Year),max(Year),by= 2),las=2) 
points(Year,(kusko.data[,'In.river']/u),pch=16,col="azure4")
arrows(Year,y0=((kusko.data[,'In.river'])+2*kusko.data[,'In.river.sd'])/u,y1=((kusko.data[,'In.river'])-2*kusko.data[,'In.river.sd'])/u,code=0,
       col="azure4",lwd=2)
arrows(Year,y0=UCIEsc/u,y1=LCIEsc[1:ny]/u,code=0, col=1,lty=5)
lines(Year,Esc/u,type='o',ljoin=1,pch=21,col=1,bg='white',lty=5)
legend('topright',legend=c('Run','Escapement','In River Run Estimate'), 
       lty=c(1,5,1),pch=c(16,21,16),col=c(1,1,"azure4"), bg=
         'white', bty ='n')

dev.off()
###############################################################################
#Figure - CV, number of projects 
##############################################################################
jpeg("CV vs. # of Projects.jpeg", width = 7, height = 4, units = "in", res = 150, pointsize = 10)

par(mar = c(5, 5, 3, 5))
projects<-rowSums(!is.na(kusko.data[,prj.names]))#number of projects operated by year. 

bp<-barplot(CV, yaxt="n", xlab="Year", ylab="Coefficient of Variation (CV)",
            axis.lty=1, ylim=c(0,0.30), names.arg=Year, las=2,xaxs = "i",xlim = c(-1, ny+10.5), cex.axis = 1.4, cex.names=.9)
box()
mtext("Number of Assessment Projects", side = 4, line = 3)

lines(x=bp, y= projects/100, lty=3, lwd=2)
lines(x=bp, y= rep(mean(CV), ny), lty=1, lwd=2) #average cv of dataset

axis(2, at=seq(0,.25,by=.05), labels=paste(100*seq(0,.25,by=.05), "%") )

axis(4, at=seq(0,.25,by=.05), labels=paste(100*seq(0,.25,by=.05)))

dev.off()
###############################################################################
#Figure Scatter Graph : Weir/Aerial vs. predicted escapement #####################################
###############################################################################
jpeg("OBS vs PREDIC by year.jpeg", width = 7, height = 7, units = "in", res = 150)

par(mfrow=c(5,5),mar = c(1.75,1.5,1.5,1.75),oma = c(3,3,3,3)) 
x <-kusko.data$In.river  #  Observed total run  
y <-(Run)  #  Model predicted 
maxx <- max(max(x,na.rm=TRUE),max(y,na.rm=TRUE))
plot(Year,y,ylab='',type='l', main='In River',ylim=c(0,maxx), cex.main =1   )
points(Year,x,pch=19,col='black')


for (i in 1:length(esc.name)){
  x <-kusko.data[,prj.names[i]]
  y <-Esc/exp(scl[i])
  maxx <- max(max(x,na.rm=TRUE),max(y,na.rm=TRUE))
  plot(Year,y,ylab='',type='l', main=esc.name[i],ylim=c(0,maxx), cex.main=1)
  points(Year,x,pch=19,col='black')
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
x <-kusko.data$In.river  #  Observed total run  
y <-(Run)  #  Model predicted 
maxx <- max(x,na.rm=TRUE)
z <-cbind(x,y)
plot(z,ylab='',col='black', main='In River',xlim=c(0,maxx),ylim=c(0,maxx),cex.main=1)
abline(0,1)

for (i in 1:length(esc.name)){
  x <-kusko.data[,prj.names[i]]
  y <-Esc/exp(scl[i])
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


par(mfrow=c(1,1),mar = c(4, 4, 1, 1))
for (i in 1:length(esc.name)){
  x <-(kusko.data[,prj.names[i]]*exp(scl[i]))/u
  z <- cbind(Year,x)
  if(i==1) plot(z[c((ny-2):ny),],pch=19,ylim=c(0,400), xlim=c(2020.5, 2023.5),xaxt="n", xlab="Year", col="azure4",
                ylab='Model Predicted Escapement (x 1,000)')
  points(z[c((ny-2):ny),], pch = 19, cex=1,col="azure4")
  
  #lines(z[c((ny-2):ny),])
  text(z[c((ny-2):ny),], labels=esc.name[i], cex= 0.5, pos=4)
}
points(Year[c((ny-2):ny)],Esc[c((ny-2):ny)]/u,col = 'black', pch = 19, 
       cex=1.5)  
arrows(Year[c((ny-2):ny)],y0=UCIEsc[c((ny-2):ny)]/1000,
       y1=LCIEsc[c((ny-2):ny)]/1000,code=0,col='black',lty=1,lwd=5)
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






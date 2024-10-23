#================================================================
# 1. Set up libraries and source functions
#===============================================================
library(R2admb)
library(RColorBrewer)
#source("KuskokwimFunctions.R")
#================================================================
# 2. Read in ADFG Monitoring data
#================================================================
setwd("C:/Users/sdlarson/OneDrive - State of Alaska/Desktop/ADMB")
kusko.data <- read.csv("KingData.csv", h=T)   ####  CSV file that was used for R Model. 
#================================================================
#================================================================
# 3. Arrange data for model input
#================================================================
#================================================================
#########################################################################################
#  3.1 Testfishery: Estimate run proporion of 1976-1983                                #   
#########################################################################################
####Make sure the RR has been run WITH headwaters before doing sensitivity ####################################################

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
init <- c(rep(log(250000),ny),rep(5,6),rep(4,14),rep(-10,3),rep(2,6),rep(2,14))	
# Lower bounds
lb <-  c(log(minrun),rep(2,6), rep(3,14),rep(-14,3),rep(-3,6),rep(-3,14))
# Upper bounds 
ub <-  c(rep(log(500000),ny),rep(7,6),rep(8,14),rep(-5,3),rep(5,6),rep(5,14))

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
#======================================================================
# Vary the level of escapmeent data to make sure the model is scaling correctly. 
#======================================================================
per<-c(0.5,0.6,0.7,0.8,0.9,1,1.1,1.2,1.3,1.4,1.5)
nsim<-length(per)


storagemat.admb<-matrix(NA,nrow=74,ncol=nsim)   #add one
fit<-rep(0,nsim)
fit2<-fit
for(i in 1:nsim)
{
  
  adfgdat2<-adfgdat
  adfgdat2$esc_w[48,]<-round(adfgdat2$esc_w[48,]*per[i],0)   #add one
  adfgdat2$esc_a[48,]<-round(adfgdat2$esc_a[48,]*per[i],0)   #add one
  
  start.parms<-list(
    log_trun = rep(12.5,ny2),
    log_wesc = rep(5,6),
    log_aesc = rep(4,14),
    log_q = rep(-11,3),
    log_cvw = 1,
    log_cva = 1,
    log_cvq = 1 
    
  )
  
  fn<-"ADFG_Rmodel_toADMB_pooled_ln"
  aout<-do_admb(fn=fn,data=adfgdat2,params=start.parms,mcmc=F,
                verbose=T,run.opts = run.control(clean_files = "none",read_files = F))
  
  adparmsout<-c(list(fn = fn),read_pars(fn))
  ad.in<-adfgdat2
  clean_admb(fn)
  storagemat.admb[,i]<-c(adparmsout$coefficients[1:48],adparmsout$coefficients[69:70],   #add one
                         adparmsout$coefficients[49:54],adparmsout$coefficients[55:68],  #add one
                         adparmsout$coefficients[71],adparmsout$coefficients[72],        #add one
                         adparmsout$loglik,adparmsout$maxgrad)
  fit[i]<-sum(is.na(adparmsout$vcov))==0
  

}
storagemat.good<-exp(storagemat.admb[1:48,fit==1])  #add one
abund<-storagemat.good[48,]   #add one
plot(abund~per, ylab="Total Run", xlab="Change (1= 0% change)", main="Total Run vs Increase/Decrease Esc")

#======================================================================
# Vary the level of harvest data to make sure the model is scaling correctly. 
#======================================================================
per<-c(0.5,0.6,0.7,0.8,0.9,1,1.1,1.2,1.3,1.4,1.5)
nsim<-length(per)


storagemat.admb<-matrix(NA,nrow=74,ncol=nsim)  #add one
fit<-rep(0,nsim)
fit2<-fit
for(i in 1:nsim)
{
  
  adfgdat2<-adfgdat
  adfgdat2$tcatch[48]<-round(adfgdat2$tcatch[48]*per[i],0)   #add one
 
  start.parms<-list(
    log_trun = rep(12.5,ny2),
    log_wesc = rep(5,6),
    log_aesc = rep(4,14),
    log_q = rep(-11,3),
    log_cvw = 1,
    log_cva = 1,
    log_cvq = 1 
    
  )
  
  fn<-"ADFG_Rmodel_toADMB_pooled_ln"
  aout<-do_admb(fn=fn,data=adfgdat2,params=start.parms,mcmc=F,
                verbose=T,run.opts = run.control(clean_files = "none",read_files = F))
  
  adparmsout<-c(list(fn = fn),read_pars(fn))
  ad.in<-adfgdat2
  clean_admb(fn)
  storagemat.admb[,i]<-c(adparmsout$coefficients[1:48],adparmsout$coefficients[69:70],   #add one
                         adparmsout$coefficients[49:54],adparmsout$coefficients[55:68],  #add one
                         adparmsout$coefficients[71],adparmsout$coefficients[72],        #add one
                         adparmsout$loglik,adparmsout$maxgrad)
  fit[i]<-sum(is.na(adparmsout$vcov))==0
  
  
}
storagemat.good<-exp(storagemat.admb[1:48,fit==1])  #add one
abund<-storagemat.good[48,]   #add one
plot(abund~per, ylab="Total Run", xlab="Change (1= 0% change)", main="Total Run vs Increase/Decrease Harvest")



########################################################################################################################################################################################################################
####################         Starting Value Evaluation                     ######################################
##################################################################################################################



######Total Run Starting Values #################################################################################################


nsim=100

startvals<-seq(100000,400000,length.out = nsim)

storagemat.admb<-matrix(NA,nrow=74,ncol=nsim)  #add one (72 in 2022?)
fit<-rep(0,nsim)
fit2<-fit
for(i in 1:nsim)
{
  
  start.parms<-list(
    log_trun = rep(log(startvals[i]),48), ##### Update number Anually #########
    log_wesc = rep(5,6),
    log_aesc = rep(4,14),
    log_q = rep(-11,3), #manipulate for commerical catch. 
    log_cvw = 1,
    log_cva = 1,
    log_cvq = 1    
  )
  
  
  fn<-"ADFG_Rmodel_toADMB_pooled_ln"
  aout<-do_admb(fn=fn,data=adfgdat,params=start.parms,mcmc=F,
                verbose=T,run.opts = run.control(clean_files = "none",read_files = F))
  
  adparmsout<-c(list(fn = fn),read_pars(fn))
  ad.in<-adfgdat
  clean_admb(fn)
  storagemat.admb[,i]<-c(adparmsout$coefficients[1:48],adparmsout$coefficients[69:70],   #add one
                         adparmsout$coefficients[49:54],adparmsout$coefficients[55:68],  #add one
                         adparmsout$coefficients[71],adparmsout$coefficients[72],        #add one
                         adparmsout$loglik,adparmsout$maxgrad)
  fit[i]<-sum(is.na(adparmsout$vcov))==0
  
}

storagemat.good<-exp(storagemat.admb[1:48,fit==1])  #add one here and below plot
ranges<-apply(storagemat.good,MARGIN=1,range)

plot(NA,xlim=c(0,48),ylim=c(0,550000),xaxt="n",xlab="year",ylab="Run Size",bty="l",xaxs="i",yaxs="i", main = "Total Run")
legend("topright",bty="n",legend=c(paste("Max Difference = ",signif(max(ranges[2,]-ranges[1,]),4),sep=""),
                                   paste("Mean Difference = ",signif(mean(ranges[2,]-ranges[1,]),4),sep="")))
for(i in 2:nrow(storagemat.good))
{
  par(new=T)
  plot(storagemat.good[,i],xaxt="n",yaxt='n',xlab="",ylab="",type='l')
  
}
temp <- storagemat.admb[74,]*fit #add one
temp[temp==0] <- NA 
hist(temp,col="grey",border="grey",
     xlab="NLL",main="")

plot(startvals,temp,
     ylab="NLL",xlab="Starting Value",pch=16)

abund<-storagemat.good[48,]   #### add one

###################Weir Starting Values#############################################################################################################

Wstartvals<-seq(0.01,10,length.out = nsim)

storagemat.admb<-matrix(NA,nrow=74,ncol=nsim)
fit<-rep(0,nsim)
fit2<-fit
for(i in 1:nsim)
{
  
  start.parms<-list(
    log_trun = rep(12.5,48), ##### Update number Anually #########
    log_wesc = rep(Wstartvals[i],6),
    log_aesc = rep(4,14),
    log_q = rep(-11,3), #manipulate for commerical catch. 
    log_cvw = 1,
    log_cva = 1,
    log_cvq = 1    
  )
  fn<-"ADFG_Rmodel_toADMB_pooled_ln"
  aout<-do_admb(fn=fn,data=adfgdat,params=start.parms,mcmc=F,
                verbose=T,run.opts = run.control(clean_files = "none",read_files = F))
  
  adparmsout<-c(list(fn = fn),read_pars(fn))
  ad.in<-adfgdat
  clean_admb(fn)
  storagemat.admb[,i]<-c(adparmsout$coefficients[1:48],adparmsout$coefficients[69:70],   #add one
                         adparmsout$coefficients[49:54],adparmsout$coefficients[55:68],  #add one
                         adparmsout$coefficients[71],adparmsout$coefficients[72],        #add one
                         adparmsout$loglik,adparmsout$maxgrad)
  fit[i]<-sum(is.na(adparmsout$vcov))==0
  
}

storagemat.good<-exp(storagemat.admb[1:48,fit==1])  #add one and below plot
ranges<-apply(storagemat.good,MARGIN=1,range)

plot(NA,xlim=c(0,48),ylim=c(0,550000),xaxt="n",xlab="year",ylab="Run Size",bty="l",xaxs="i",yaxs="i", main = "Weir")
legend("topright",bty="n",legend=c(paste("Max Difference = ",signif(max(ranges[2,]-ranges[1,]),4),sep=""),
                                   paste("Mean Difference = ",signif(mean(ranges[2,]-ranges[1,]),4),sep="")))
for(i in 2:nrow(storagemat.good))
{
  par(new=T)
  plot(storagemat.good[,i],xaxt="n",yaxt='n',xlab="",ylab="",type='l')
  
}
temp <- storagemat.admb[74,]*fit    # update
temp[temp==0] <- NA 
hist(temp,col="grey",border="grey",
     xlab="NLL",main="")

plot(Wstartvals,temp,
     ylab="NLL",xlab="Starting Value",pch=16)

weirs<-storagemat.good[48,]  #add one

###################Aerial Starting Values#############################################################################################################

Astartvals<-seq(0.01,10,length.out = nsim)

storagemat.admb<-matrix(NA,nrow=74,ncol=nsim)   #add one
fit<-rep(0,nsim)
fit2<-fit
for(i in 1:nsim)
{
  
  start.parms<-list(
    log_trun = rep(12.5,48), ##### Update number Anually #########
    log_wesc = rep(5,6),
    log_aesc = rep(Astartvals[i],14),
    log_q = rep(-11,3), #manipulate for commerical catch. 
    log_cvw = 1,
    log_cva = 1,
    log_cvq = 1      
  )
  fn<-"ADFG_Rmodel_toADMB_pooled_ln"
  aout<-do_admb(fn=fn,data=adfgdat,params=start.parms,mcmc=F,
                verbose=T,run.opts = run.control(clean_files = "none",read_files = F))
  
  adparmsout<-c(list(fn = fn),read_pars(fn))
  ad.in<-adfgdat
  clean_admb(fn)
  storagemat.admb[,i]<-c(adparmsout$coefficients[1:48],adparmsout$coefficients[69:70],   #add one
                         adparmsout$coefficients[49:54],adparmsout$coefficients[55:68],  #add one
                         adparmsout$coefficients[71],adparmsout$coefficients[72],        #add one
                         adparmsout$loglik,adparmsout$maxgrad)
  fit[i]<-sum(is.na(adparmsout$vcov))==0
  
}

storagemat.good<-exp(storagemat.admb[1:48,fit==1]) #add one here and below
ranges<-apply(storagemat.good,MARGIN=1,range)

plot(NA,xlim=c(0,48),ylim=c(0,550000),xaxt="n",xlab="year",ylab="Run Size",bty="l",xaxs="i",yaxs="i",main = "Air")
legend("topright",bty="n",legend=c(paste("Max Difference = ",signif(max(ranges[2,]-ranges[1,]),4),sep=""),
                                   paste("Mean Difference = ",signif(mean(ranges[2,]-ranges[1,]),4),sep="")))
for(i in 2:nrow(storagemat.good))
{
  par(new=T)
  plot(storagemat.good[,i],xaxt="n",yaxt='n',xlab="",ylab="",type='l')
  
}

temp <- storagemat.admb[74,]*fit #update
temp[temp==0] <- NA 
hist(temp,col="grey",border="grey",
     xlab="NLL",main="")

plot(Astartvals,temp,
     ylab="NLL",xlab="Starting Value",pch=16)
aerial<-storagemat.good[48,]  #add one


###################comm Starting Values#############################################################################################################

commstartvals<-seq(-20,1,length.out = nsim)

storagemat.admb<-matrix(NA,nrow=74,ncol=nsim)   #add one
fit<-rep(0,nsim)
fit2<-fit
for(i in 1:nsim)
{
  
  start.parms<-list(
    log_trun = rep(12.5,48), ##### Update number Anually #########
    log_wesc = rep(5,6),
    log_aesc = rep(4,14),
    log_q = rep(commstartvals[i],3), #manipulate for commerical catch. 
    log_cvw = 1,
    log_cva = 1,
    log_cvq = 1      
  )
  fn<-"ADFG_Rmodel_toADMB_pooled_ln"
  aout<-do_admb(fn=fn,data=adfgdat,params=start.parms,mcmc=F,
                verbose=T,run.opts = run.control(clean_files = "none",read_files = F))
  
  adparmsout<-c(list(fn = fn),read_pars(fn))
  ad.in<-adfgdat
  clean_admb(fn)
  storagemat.admb[,i]<-c(adparmsout$coefficients[1:48],adparmsout$coefficients[69:70],   #add one
                         adparmsout$coefficients[49:54],adparmsout$coefficients[55:68],  #add one
                         adparmsout$coefficients[71],adparmsout$coefficients[72],        #add one
                         adparmsout$loglik,adparmsout$maxgrad)
  fit[i]<-sum(is.na(adparmsout$vcov))==0
}

storagemat.good<-exp(storagemat.admb[1:48,fit==1]) #add one and below
ranges<-apply(storagemat.good,MARGIN=1,range)

plot(NA,xlim=c(0,48),ylim=c(0,550000),xaxt="n",xlab="year",ylab="Run Size",bty="l",xaxs="i",yaxs="i",main = "Comm")
legend("topright",bty="n",legend=c(paste("Max Difference = ",signif(max(ranges[2,]-ranges[1,]),4),sep=""),
                                   paste("Mean Difference = ",signif(mean(ranges[2,]-ranges[1,]),4),sep="")))
for(i in 2:nrow(storagemat.good))
{
  par(new=T)
  plot(storagemat.good[,i],xaxt="n",yaxt='n',xlab="",ylab="",type='l')
  
}
temp <- storagemat.admb[74,]*fit   #update
temp[temp==0] <- NA 
hist(temp,col="grey",border="grey",
     xlab="NLL",main="")

plot(commstartvals,temp,
     ylab="NLL",xlab="Starting Value",pch=16)

comm<-storagemat.good[48,]   #add one



###################cv Starting Values#############################################################################################################
cvstartvals<-seq(-20,20,length.out = nsim)

storagemat.admb<-matrix(NA,nrow=74,ncol=nsim)  #add one
fit<-rep(0,nsim)
fit2<-fit
for(i in 1:nsim)
{
  
  start.parms<-list(
    log_trun = rep(12.5,48), ##### Update number Anually #########
    log_wesc = rep(5,6),
    log_aesc = rep(4,14),
    log_q = rep(-11,3), #manipulate for commerical catch. 
    log_cvw = cvstartvals[i],
    log_cva = cvstartvals[i],
    log_cvq = cvstartvals[i]    
  )
  fn<-"ADFG_Rmodel_toADMB_pooled_ln"
  aout<-do_admb(fn=fn,data=adfgdat,params=start.parms,mcmc=F,
                verbose=T,run.opts = run.control(clean_files = "none",read_files = F))
  
  adparmsout<-c(list(fn = fn),read_pars(fn))
  ad.in<-adfgdat
  clean_admb(fn)
  storagemat.admb[,i]<-c(adparmsout$coefficients[1:48],adparmsout$coefficients[69:70],   #add one
                         adparmsout$coefficients[49:54],adparmsout$coefficients[55:68],  #add one
                         adparmsout$coefficients[71],adparmsout$coefficients[72],        #add one
                         adparmsout$loglik,adparmsout$maxgrad)
  fit[i]<-sum(is.na(adparmsout$vcov))==0
  
}

storagemat.good<-exp(storagemat.admb[1:48,fit==1]) #add one here and below
ranges<-apply(storagemat.good,MARGIN=1,range)

plot(NA,xlim=c(0,48),ylim=c(0,550000),xaxt="n",xlab="year",ylab="Run Size",bty="l",xaxs="i",yaxs="i",main = "CV")
legend("topright",bty="n",legend=c(paste("Max Difference = ",signif(max(ranges[2,]-ranges[1,]),4),sep=""),
                                   paste("Mean Difference = ",signif(mean(ranges[2,]-ranges[1,]),4),sep="")))

for(i in 2:nrow(storagemat.good))
{
  par(new=T)
  plot(storagemat.good[,i],xaxt="n",yaxt='n',xlab="",ylab="",type='l')
  
}

temp <- storagemat.admb[74,]*fit  #add one
temp[temp==0] <- NA 
hist(temp,col="grey",border="grey",
     xlab="NLL",main="")

plot(cvstartvals,temp,
     ylab="NLL",xlab="Starting Value",pch=16)

cv<-storagemat.good[48,]  #add one
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                             
####################################################################################################################################################


###################### Project Removal ############################

###    Weirs    #####

storagemat.admb<-matrix(NA,nrow=74,ncol=length(adfgdat$esc_w))  #add one
fit<-rep(0,length(adfgdat$esc_w))
fit2<-fit
for(i in 1:length(adfgdat$esc_w))
{
  
  adfgdat3<-adfgdat
  adfgdat3$esc_w[48,i]<-0  #add one
 
  start.parms<-list(
    log_trun = rep(12.5,ny2),
    log_wesc = rep(5,6),
    log_aesc = rep(4,14),
    log_q = rep(-11,3),
    log_cvw = 1,
    log_cva = 1,
    log_cvq = 1 
    
  )
  
  fn<-"ADFG_Rmodel_toADMB_pooled_ln"
  aout<-do_admb(fn=fn,data=adfgdat3,params=start.parms,mcmc=F,
                verbose=T,run.opts = run.control(clean_files = "none",read_files = F))
  
  adparmsout<-c(list(fn = fn),read_pars(fn))
  ad.in<-adfgdat3
  clean_admb(fn)
  storagemat.admb[,i]<-c(adparmsout$coefficients[1:48],adparmsout$coefficients[69:70],   #add one
                         adparmsout$coefficients[49:54],adparmsout$coefficients[55:68],  #add one
                         adparmsout$coefficients[71],adparmsout$coefficients[72],        #add one
                         adparmsout$loglik,adparmsout$maxgrad)
  fit[i]<-sum(is.na(adparmsout$vcov))==0
  
  
}
storagemat.good<-exp(storagemat.admb[1:48,fit==1])  #add one
Single.Weir<-storagemat.good[48,]   #add one

sean<-adparmsout$se

###    Aerial    #####

storagemat.admb<-matrix(NA,nrow=74,ncol=length(adfgdat$esc_a))  #add one
fit<-rep(0,length(adfgdat$esc_a))
fit2<-fit
for(i in 1:length(adfgdat$esc_a))
{
  
  adfgdat3<-adfgdat
  adfgdat3$esc_a[48,i]<-0  #add one
  
  start.parms<-list(
    log_trun = rep(12.5,ny2),
    log_wesc = rep(5,6),
    log_aesc = rep(4,14),
    log_q = rep(-11,3),
    log_cvw = 1,
    log_cva = 1,
    log_cvq = 1 
    
  )
  
  fn<-"ADFG_Rmodel_toADMB_pooled_ln"
  aout<-do_admb(fn=fn,data=adfgdat3,params=start.parms,mcmc=F,
                verbose=T,run.opts = run.control(clean_files = "none",read_files = F))
  
  adparmsout<-c(list(fn = fn),read_pars(fn))
  ad.in<-adfgdat3
  clean_admb(fn)
  storagemat.admb[,i]<-c(adparmsout$coefficients[1:48],adparmsout$coefficients[69:70],   #add one
                         adparmsout$coefficients[49:54],adparmsout$coefficients[55:68],  #add one
                         adparmsout$coefficients[71],adparmsout$coefficients[72],        #add one
                         adparmsout$loglik,adparmsout$maxgrad)
  fit[i]<-sum(is.na(adparmsout$vcov))==0
  
  
}
storagemat.good<-exp(storagemat.admb[1:48,fit==1])  #add one
Single.Air<-storagemat.good[48,]  #add one


###    No Weirs  ##### works in 2021
storagemat.admb<-matrix(NA,nrow=74,ncol=1)   #add one
fit<-rep(0,1)
fit2<-fit

  
  adfgdat3<-adfgdat
  adfgdat3$esc_w[48,]<-0   #add one
  
  start.parms<-list(
    log_trun = rep(12.5,ny2),
    log_wesc = rep(5,6),
    log_aesc = rep(4,14),
    log_q = rep(-11,3),
    log_cvw = 1,
    log_cva = 1,
    log_cvq = 1 
    
  )
  
  fn<-"ADFG_Rmodel_toADMB_pooled_ln"
  aout<-do_admb(fn=fn,data=adfgdat3,params=start.parms,mcmc=F,
                verbose=T,run.opts = run.control(clean_files = "none",read_files = F))
  
  adparmsout<-c(list(fn = fn),read_pars(fn))
  ad.in<-adfgdat3
  clean_admb(fn)
  storagemat.admb[,1]<-c(adparmsout$coefficients[1:48],adparmsout$coefficients[69:70],   #add one
                         adparmsout$coefficients[49:54],adparmsout$coefficients[55:68],  #add one
                         adparmsout$coefficients[71],adparmsout$coefficients[72],        #add one
                         adparmsout$loglik,adparmsout$maxgrad)
  fit[i]<-sum(is.na(adparmsout$vcov))==0
  

storagemat.good<-exp(storagemat.admb[1:48])  #add one  ,fit==1
No.Weirs<-storagemat.good[[48]] #aerial only for Figure 6                             ################################################### FIGURE
adparmsout$se  #SE for table 6

###    No Air  ##### works in 2021
storagemat.admb<-matrix(NA,nrow=74,ncol=1)  #add one
fit<-rep(0,1)
fit2<-fit


adfgdat3<-adfgdat
adfgdat3$esc_a[48,]<-0   #add one

start.parms<-list(
  log_trun = rep(12.5,ny2),
  log_wesc = rep(5,6),
  log_aesc = rep(4,14),
  log_q = rep(-11,3),
  log_cvw = 1,
  log_cva = 1,
  log_cvq = 1 
  
)

fn<-"ADFG_Rmodel_toADMB_pooled_ln"
aout<-do_admb(fn=fn,data=adfgdat3,params=start.parms,mcmc=F,
              verbose=T,run.opts = run.control(clean_files = "none",read_files = F))

adparmsout<-c(list(fn = fn),read_pars(fn))
ad.in<-adfgdat3
clean_admb(fn)
storagemat.admb[,1]<-c(adparmsout$coefficients[1:48],adparmsout$coefficients[69:70],   #add one
                       adparmsout$coefficients[49:54],adparmsout$coefficients[55:68],  #add one
                       adparmsout$coefficients[71],adparmsout$coefficients[72],        #add one
                       adparmsout$loglik,adparmsout$maxgrad)
fit[i]<-sum(is.na(adparmsout$vcov))==0



storagemat.good<-exp(storagemat.admb[1:48])       #add one    ,fit==1
No.Air<-storagemat.good[[48]] #weir only for Figure 6    #add one
adparmsout$se  #SE for table 6

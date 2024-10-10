#########################################################################################
#   Kuskokwim Chinook Salmon Run Reconstruction R-code
#
#   Written by  H. Hamaazaki 
#   Date          
#   
#########################################################################################

#########################################################################################

#########################################################################################
#  1.0  Initialize working Environment                                                  #   
#########################################################################################
rm(list=ls(all=TRUE))
# Enter the name of data file  
data_file <- 'Kusko_RR_Input_March_2_17.csv'

kusko.data <- read.csv(data_file,header=T, na.string='')

#########################################################################################
#  2.2  Testfishery: Estimate run proporion of 1976-1983                                #   
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

#########################################################################################
#  2.3  Rearange fishing effort and harvest data catch 0 to NA                          #   
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
ccat[ccat == 0] <- NA
# Extract weekly commercial est data
creg <-kusko.data[substr(names(kusko.data),1,3)=='cfw']
# combine week 8, 9  and drop 
creg[,6] <- pmax(creg[,6],creg[,7])
creg <- creg[,-7]

#########################################################################################
#  2.4  Recalculate Inriver data                           #   
#########################################################################################
# Extract Inriver data 
inr <-kusko.data[substr(names(kusko.data),1,3)=='In.']
# Calculate CV
inr$cv <- inr$In.river.sd/inr$In.river

#########################################################################################
#  2.5  Calculate Others                                              #   
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
#  2.4 Construct dataset used for likelihood modeling                                               #   
#########################################################################################
kusko.like.data <- as.matrix(cbind(tcatch,inr,esc,testf[3:8],ccat,ceff,creg))

nb.likelihood <- function(theta,likedat,ny){
### Total run  ###################################################################
  	totrun <- exp(theta[1:ny])
### Weir slope parameters ###################################################################
  	w.kwe <- exp(theta[ny+1])
  	w.tul <- exp(theta[ny+2])
  	w.geo <- exp(theta[ny+3])
  	w.kog <- exp(theta[ny+4])
  	w.tat <- exp(theta[ny+5])
  	w.tak <- exp(theta[ny+6])
### Aerial slope parameters ###################################################################  	
  	a.kwe <- exp(theta[ny+7])
  	a.kis <- exp(theta[ny+8])
  	a.tul <- exp(theta[ny+9])
    a.sla <- exp(theta[ny+10])	
  	a.kip <- exp(theta[ny+11])
  	a.ank <- exp(theta[ny+12])
  	a.hlk <- exp(theta[ny+13])
    a.osk <- exp(theta[ny+14])
  	a.hlt <- exp(theta[ny+15])
  	a.che <- exp(theta[ny+16])
  	a.gag <- exp(theta[ny+17])
  	a.pit <- exp(theta[ny+18])
  	a.ber <- exp(theta[ny+19])
    a.slp <- exp(theta[ny+20])

### Catch coefficient parameters ###################################################################
# catchability coefficient Unrestricted
  	q1 <- exp(theta[ny+21])
# catchability coefficient Restricted
  	q2 <- exp(theta[ny+22])
# catchability coefficient Center Core monofilament	
  	q3 <- exp(theta[ny+23])

### Overdispersion Parameters ########################################################################
## Weir #########################################################################################  	
	r.kwe <- exp(theta[ny+24])
  	r.tul <- exp(theta[ny+25])
  	r.geo <- exp(theta[ny+26])
  	r.kog <- exp(theta[ny+27])
  	r.tat <- exp(theta[ny+28])
  	r.tak <- exp(theta[ny+29])
### Aerial ###################################################################  	
  	ra.kwe <- exp(theta[ny+30])
  	ra.kis <- exp(theta[ny+31])
  	ra.tul <- exp(theta[ny+32])
    ra.sla <- exp(theta[ny+33])	
  	ra.kip <- exp(theta[ny+34])
  	ra.ank <- exp(theta[ny+35])
  	ra.hlk <- exp(theta[ny+36])
    ra.osk <- exp(theta[ny+37])
  	ra.hlt <- exp(theta[ny+38])
  	ra.che <- exp(theta[ny+39])
  	ra.gag <- exp(theta[ny+40])
  	ra.pit <- exp(theta[ny+41])
  	ra.ber <- exp(theta[ny+42])
    ra.slp <- exp(theta[ny+43])	
### Likelihood model ###################################################################  	     	
    tfw  <-  rep(0,6)
    tfa  <-  rep(0,14)
    tft  <-  0
    tfc  <-  0
	esc  <-  totrun-likedat[,1]
	
#### Definie the negative binomial function ############################################
nblike <- function(obs,r,est){
	lgamma(obs+r)-lgamma(obs+1)-lgamma(r)+r*log(r/(est+r))+obs*log(est/(est+r))
	}
#### Weir likelifhood ###################################################################    
    tfw[1]  <-  -sum(nblike(likedat[,5],r.kwe,esc/w.kwe),na.rm=T)  
    tfw[2]  <-  -sum(nblike(likedat[,6],r.tul,esc/w.tul),na.rm=T)   
    tfw[3]  <-  -sum(nblike(likedat[,7],r.geo,esc/w.geo),na.rm=T)
    tfw[4]  <-  -sum(nblike(likedat[,8],r.kog,esc/w.kog),na.rm=T)
    tfw[5]  <-  -sum(nblike(likedat[,9],r.tat,esc/w.tat),na.rm=T)
    tfw[6]  <-  -sum(nblike(likedat[,10],r.tak,esc/w.tak),na.rm=T)
#### Aerial likelihood ####################################################################    
    tfa[1]  <-  -sum(nblike(likedat[,11],ra.kwe,esc/a.kwe),na.rm=T)
    tfa[2]  <-  -sum(nblike(likedat[,12],ra.kis,esc/a.kis),na.rm=T)
    tfa[3]  <-  -sum(nblike(likedat[,13],ra.tul,esc/a.tul),na.rm=T)
    tfa[4]  <-  -sum(nblike(likedat[,14],ra.sla,esc/a.sla),na.rm=T)
    tfa[5]  <-  -sum(nblike(likedat[,15],ra.kip,esc/a.kip),na.rm=T)
    tfa[6]  <-  -sum(nblike(likedat[,16],ra.ank,esc/a.ank),na.rm=T)
    tfa[7]  <-  -sum(nblike(likedat[,17],ra.hlk,esc/a.hlk),na.rm=T)
    tfa[8]  <-  -sum(nblike(likedat[,18],ra.osk,esc/a.osk),na.rm=T)
    tfa[9]  <-  -sum(nblike(likedat[,19],ra.hlt,esc/a.hlt),na.rm=T)
    tfa[10]  <-  -sum(nblike(likedat[,20],ra.che,esc/a.che),na.rm=T)
    tfa[11]  <-  -sum(nblike(likedat[,21],ra.gag,esc/a.gag),na.rm=T)
    tfa[12]  <-  -sum(nblike(likedat[,22],ra.pit,esc/a.pit),na.rm=T)
    tfa[13]  <-  -sum(nblike(likedat[,23],ra.ber,esc/a.ber),na.rm=T)
    tfa[14]  <-  -sum(nblike(likedat[,24],ra.slp,esc/a.slp),na.rm=T)
    
#### Inriver  Normal likelifhood #######################################################################
    tft  <-  0.5*sum((likedat[,2]-totrun)^2/(likedat[,3])^2,na.rm=T)

#### Weekly Catch likelifhood ###################################################################  
# Calculated estimated run by week 
wk.est <- likedat[,25:30]*totrun
#################################################################################################
#   Calculate likelihoood for unrestricted 
#################################################################################################
# Extract all mesh regulation year/week 
unr <- likedat[,43:48]
# Keep unrestricted mesh regulation year/week 1: indicate unrestricted period
unr[unr != 1] <- NA
# Observed Effort
# Keep only Effort of Unrestricted 
unr.eff <- likedat[,37:42]*unr
# Rmove all NA 
unr.eff <- unr.eff[!is.na(unr.eff)]

# Observed harvest
# Keep only Effort of Unrestricted 
unr.h <- likedat[,31:36]*unr
# Rmove all NA 
unr.h <- unr.h[!is.na(unr.h)]

# Estimated 
# Keep only Effort of Unrestricted 
unr.wk <- wk.est*unr
# Rmove all NA 
unr.wk <- unr.wk[!is.na(unr.wk)]

# likelihood for Unrestricted
   tf1 <- 0.5*length(unr.eff)*log(sum((log(unr.eff)-log(-log(1-unr.h/unr.wk)/q1))^2,na.rm=T))
  
#################################################################################################
#   Calculate likelihoood for restricted 
#################################################################################################
# Extract restricted mesh period 
# Extract all mesh regulation year/week 
r <- likedat[,43:48]
# Keep unrestricted mesh regulation year/week 2: indicate restricted periods
r[r != 2] <- NA
# Change it to 1 
r[r == 2] <- 1
# Observed effort
# Keep only Effort of Restricted 
r.eff <- likedat[,37:42]*r
# Rmove all NA 
r.eff <- r.eff[!is.na(r.eff)]

# Observed harvest
# Keep only Effort of Restricted 
r.h <- likedat[,31:36]*r
# Rmove all NA 
r.h <- r.h[!is.na(r.h)]

# Estimated
# Keep only Effort of Unrestricted 
r.wk <- wk.est*r
# Rmove all NA 
r.wk <- r.wk[!is.na(r.wk)]
 
# likelihood for Unrestricted
   tf2 <- 0.5*length(r.eff)*log(sum((log(r.eff)-log(-log(1-r.h/r.wk)/q2))^2,na.rm=T))

#################################################################################################
#   Calculate likelihoood for Monofilament
#################################################################################################
# Extract Monfilament periods
# Extract all mesh regulation year/week (This is taking only 3-6 weeks
m <- likedat[,43:48]
# Keep monofilament mesh regulation year/week 3: indicate monofilament peiriods
m[(m != 3)&(m != 5)] <- NA
# Change it to 1 
m[!is.na(m)] <- 1
# Observed effort
# Keep only Effort of Restricted 
m.eff <- likedat[,37:42]*m
# Rmove all NA 
m.eff <- m.eff[!is.na(m.eff)]

# Observed harvest
# Keep only Effort of Restricted 
m.h <- likedat[,31:36]*m
# Rmove all NA 
m.h <- m.h[!is.na(m.h)]

# Estimated
# Keep only Effort of Restricted
m.wk <- wk.est*m
# Rmove all NA 
m.wk <- m.wk[!is.na(m.wk)]

   tf3 <- 0.5*length(m.eff)*log(sum((log(m.eff)-log(-log(1-ifelse(m.h/m.wk<1,m.h/m.wk,0.999))/q3))^2,na.rm=T))

tfc <-sum(tf1,tf2,tf3)

#### Likelihood calculation  ###################################################################    
loglink  <- sum(sum(tfw),sum(tfa),tft,tfc,na.rm=T)
return(loglink) 
}

#########################################################################################
#  3.1 Set Initial value and boundaries                                                 #   
#########################################################################################
# Initial starting point
init <- c(rep(log(250000),ny),rep(5,6),rep(4,14),rep(-10,3),rep(2,6),rep(2,14))	
# Lower bounds
lb <-  c(log(minrun),rep(2,6), rep(3,14),rep(-14,3),rep(-3,6),rep(-3,14))
# Upper bounds 
ub <-  c(rep(log(500000),ny),rep(7,6),rep(8,14),rep(-5,3),rep(5,6),rep(5,14))

#########################################################################################
#  3.3 Run likelihood model                                                      #   
#########################################################################################
ptm <- proc.time()
nll <- optim(par=init,fn=nb.likelihood,method="L-BFGS-B",lower=lb, upper = ub, control = list(maxit=1000),likedat=kusko.like.data, ny=ny, hessian = T)
min_NLL <- nll$value
proc.time() - ptm
nll$convergence
Rprof()
nll$par
nll$value

#################################################################################################
# 3.4  Calculate Wald Confidence Interval 
#################################################################################################
#1: Hessian Matrix
hessian_obs <- nll$hessian

log_est_obs <- nll$par
est_obs <- exp(log_est_obs)
# Create a variance-covariance matrix
#if (method==1) {
#    hessian_obs <- hessian_obs[1:(dim(hessian_obs)[1]),1:(dim(hessian_obs)[2])]
#    }    
var_covar_mat_obs <- solve(hessian_obs)
# Pull out diagonal nu
log_var_obs <- diag(var_covar_mat_obs)
# Calculate standard error
log_std_err_obs <- sqrt(log_var_obs)
#if (method==1) {
#    log_std_err_obs <- c(log_std_err_obs,0,0,0)
#    }
upper95CI <- exp(log_est_obs + 1.96*log_std_err_obs)
lower95CI <- exp(log_est_obs - 1.96*log_std_err_obs)
labelT <- length(ny)

for (i in 1:ny){
labelT[i] <- paste('Run',1975+i)
}
labelT <- c(labelT,names(esc),'q1','q2','q3',names(esc))
output <- data.frame(parameter=labelT,mean=exp(nll$par),lower95CI=lower95CI,lower95CI=upper95CI)

#################################################################################################
#From this point on added by Zachary Liller
#Not included in published code
#################################################################################################

# Calculate Total Run, Esc, and CI
totrun <- exp(nll$par[1:ny])
trunuci <- upper95CI[1:ny]   
trunlci <- lower95CI[1:ny]   
escp  <- totrun - tcatch
escuci <- trunuci - tcatch
esclci <- trunlci - tcatch

output2 <- data.frame(parameter=labelT[1:41],totrun=totrun,trunlci=trunlci,trunuci=trunuci,escp =escp,esclci=esclci,escuci=escuci)

write.csv(output,file='Parameter_Summary.csv')
write.csv(output2,file='Run_Esc_Summary.csv')

################################################################################
# Graphics for publication
################################################################################

windows(record=TRUE)
 
###############################################################################
#Figure - Total Run, Total Escapement, Scalars
##############################################################################
par(mfrow=c(1,1),mar = c(4, 4, 1, 1)) 
year <- seq(1976,(1976+ny-1))
#Graph
plot(year,totrun/1000,type='o',ljoin=1,ylab = 'Total run x 1000',pch=16,yaxt='n'
     , xaxt='n',ylim=c(0,550))
arrows(year,y0= trunuci/1000,y1=trunlci/1000,code=0)
axis(side=2, at = seq(0,550,by= 50),las=2) 
axis(side=1, at = seq(1976,1976+ny,by= 2), las=2) 
points(year,(inr[,1])/1000,pch=16,col='red')
arrows(year,y0=((inr[,1])+2*inr[,2])/1000,y1=((inr[,1])-2*inr[,2])/1000,code=0,
       col='red',lwd=2)
arrows(year,y0=escuci/1000,y1=esclci/1000,code=0, col=1,lty=5)
lines(year,escp/1000,type='o',ljoin=1,pch=21,yaxt='n', xaxt='n',ylim=c(0,550),
      col=1,bg='white',lty=5)
legend('topright',legend=c('Run','Escapement','In River Run estimate'), 
       lty=c(1,5,1),pch=c(16,21,16),col=c(1,1,2), bg=
         'white', bty ='n')

###############################################################################
#Figure - Scatter Graph : Weir/Aerial vs. predicted escapement
###############################################################################
par(mfrow=c(5,5),mar = c(1.75,1.5,1.5,1.75),oma = c(3,3,3,3),cex=0.6) 
x <-kusko.data[6]
y <-(totrun)
maxx <- max(x,na.rm=TRUE)
z <-cbind(x,y)
plot(z,ylab='',col='red', main='In River',xlim=c(0,maxx),ylim=c(0,maxx))
abline(0,1)
esc.name <- c('Kwethluk Weir','Tuluksak Weir','George Weir','Kogrukluk Weir',
              'Tatlawiksuk Weir','Takotna Weir','Kwethluk Aerial','Kisaralik Aerial',
              'Tuluksak Aerial','Salmon (Aniak) Aerial','Kipchuk Aerial','Aniak Aerial',
              'Holokuk Aeriai','Oskawalik Aerial','Holitna Aerial','Cheeneetnuk Aerial',
              'Gagaryah Aerial','Pitka Aerial',
              'Bear Aerial','Salmon (Pitka) Aerial')
for (i in 1:20){
  x <-kusko.data[7+i]
  y <-escp/exp(nll$par[ny+i])
  maxx <- max(x,na.rm=TRUE)
  z <-cbind(x,y)
  plot(z,ylab='',col='red', main=esc.name[i], xlim=c(0,maxx),ylim=c(0,maxx))
  points(z[ny,],col='red', pch = 19, cex=1.5)
  RMSE <- round(sqrt(sum((log(x)-log(y))^2,na.rm=T)/sum(!is.na(x))),2)
  tex <- paste('RMSE ',RMSE,sep = "")
  legend('topleft',tex,bty='n',xjust=0,inset=-0.05)
  abline(0,1)
}
#  Add Texts to Figure
mtext("Observed vs. Estimated Counts", side = 3, line = 0, outer = TRUE)
mtext('Predicted', side = 2, line = 1, outer = TRUE)
mtext("Observed", side = 1, line = 1, outer = TRUE)

###############################################################################
#Figure - Esc. by project overlayed with drainage escapement, 2014-2016
###############################################################################

par(mfrow=c(1,1),mar = c(4, 4, 1, 1))
for (i in 1:20){
  x <-(kusko.data[7+i]*exp(nll$par[ny+i]))/1000
  z <- cbind(year,x)
  if(i==1) plot(z[c((ny-2):ny),],pch=19,ylim=c(0,450), xlim=c(2013.5, 2017),xaxt="n",
                ylab='Model Predicted Escapement (x 1000)')
  points(z[c((ny-2):ny),], pch = 19, cex=1)
  #lines(z[c((ny-2):ny),])
  text(z[c((ny-2):ny),], labels=esc.name[i], cex= 0.7, pos=4)
}
points(year[c((ny-2):ny)],escp[c((ny-2):ny)]/1000,col = 'red', pch = 19, 
       cex=1.5)  
arrows(year[c((ny-2):ny)],y0=escuci[c((ny-2):ny)]/1000,
       y1=esclci[c((ny-2):ny)]/1000,code=0,col=2,lty=1,lwd=3)
axis(1,at=c(2014,2015,2016))


#'###########################################################################
#   Kuskokwim Chinook Salmon Run reconstruction data construction code
#   Kusko_Create_data.R 
#   Version 2.0 
#   Written by:   T. Hamaazaki 
#   Date: 07/03/2023 
#   This program creates an input data for Kuskokowim Run Reconstructiton model
#   The program use two data sources 
#   data_file_1:  Harvest, Run, Escapement, fishery data
#   data_file_2:  Harvest and Escapement Age Composition data        
#'###########################################################################

#'###########################################################################
#   1.0 Souce file locations: 
#'###########################################################################
library(reshape2)
# Set Base directory 

#'------------------------------------------------------------------------------
#  2.1  Read Escapement and Harvest CSV Data                                                        
#'------------------------------------------------------------------------------
setwd(Data_dir)

kusko.data <- read.csv(data_file1,header=TRUE, na.string='')

#'------------------------------------------------------------------------------
#  2.2  Import Age CSV Data                                                        
#'------------------------------------------------------------------------------
age.data <- read.csv(data_file2,header=TRUE, na.string='')
# Extract age data 
age.class <- age.data[!names(age.data) %in% c('Efn_H','Efn_E')]
# Change wide to long: This will put column name as "variable" and age comp as "value"
age.class <- melt(age.class, id.vars='Year')
# Extract Havest(H) vs Escapement (E) from the variable column 
age.class$loc <- substr(age.class$variable,1,1)
# Extract age from the variable column
age.class$at <- as.numeric(substr(age.class$variable,3,5))
# Convert European age notation to simple age 
age.class$Age <- with(age.class, floor(at)+10*(at-floor(at))+1)
age.class$Age[age.class$Age<4] <-4
age.class$Age[age.class$Age>7] <-7
# Sum proportion by age 
age.class.sum <- aggregate(value~loc+Age+Year, sum, data=age.class)
# Change back long to Wide format 
age.class  <- dcast(age.class.sum,Year~loc+Age,value.var='value')
# Get Harvest age (select column starting from H)
age.class.h <- age.class[substr(names(age.class),1,1)=='H']
# Get Standardized proprion 
age.class.h <- age.class.h/rowSums(age.class.h)  
# Extract Escapement age (select column starting from E) 
age.class.e <- age.class[substr(names(age.class),1,1)=='E']
# Get Standardized proprion 
age.class.e <- age.class.e/rowSums(age.class.e)  
# Select age sample size for escapement and harvest 
efn_H <- age.data$Efn_H
efn_E <- age.data$Efn_E

#'==============================================================================
#  2.21 Data summary Figure   
#'==============================================================================
# Read lookup table
lookup <- read.csv(data_file3,header = TRUE, na.string='')
# Extract Weir escapement data
# Mark-recapcutre data 
mr <- lookup[which(lookup$River =='MR'),]
mrs <- kusko.data[,names(kusko.data) %in% mr$CSV]
#: tmr (lower river mr), mr(Kalskag mr), amr(aniak mr), hmr (Holitna mr)
# Weir data
w.esc <- kusko.data[substr(names(kusko.data),1,2)=='w.']
# Aerial data 
a.esc <- kusko.data[substr(names(kusko.data),1,2)=='a.']
# Sonar data 
sonar <- kusko.data[,c('Sonar')]

# Extract Upriver data 
up <- lookup[which(lookup$River =='Up'),]
Esc.fig.u1 <- kusko.data[,names(kusko.data) %in% up$CSV]
Esc.fig.u2 <- kusko.data[,names(kusko.data) %in% up$CSV]
mid <- lookup[which(lookup$River =='Mid'),]
Esc.fig.m1 <- kusko.data[,names(kusko.data) %in% mid$CSV]
Esc.fig.m2 <- kusko.data[,names(kusko.data) %in% mid$CSV]
low <- lookup[which(lookup$River =='Low'),]
Esc.fig.l1 <- kusko.data[,names(kusko.data) %in% low$CSV]
Esc.fig.l2 <- kusko.data[,names(kusko.data) %in% low$CSV]
M.fig.s1 <- cbind(mrs,kusko.data$Sonar)
M.fig.s2 <- cbind(mrs,kusko.data$Sonar)
# Erase data 
Esc.fig.u2[,] <- NA
Esc.fig.m2[,] <- NA
Esc.fig.l2[,] <- NA
M.fig.s2[,-3] <-NA
year <- kusko.data$Year
fig.dat1 <- list(Esc.fig.u1,Esc.fig.m1,Esc.fig.l1,M.fig.s1)
fig.dat2 <- list(Esc.fig.u2,Esc.fig.m2,Esc.fig.l2,M.fig.s2)
n.u <- dim(Esc.fig.u1)[2]
n.m <- dim(Esc.fig.m1)[2]
n.l <- dim(Esc.fig.l1)[2]
n.s <- dim(M.fig.s1)[2]
n.data <- c(n.u, n.m, n.l, n.s)
sn.data <- sum(n.data)+4


riv <- c('Upper','Middle','Lower','Mainstem')
label.u <- up$Label
label.m <- mid$Label
label.l <- low$Label	 
label.s <- c(mr$Label,'Bethel Sonar')	
label.dat <- list(label.u,label.m,label.l,label.s)  

par(mfrow=c(1,1),mar=c(2,2,2,9),oma=c(3,1,1,1))
plot(year,Esc.fig.u1[,1],xlim=c(min(year)-1,max(year+1)),ylim=c(0,sn.data),xaxs='i',yaxt='n',
       type="n",xlab="",ylab="",main="",cex.main=1.)  
st <- 0
for (i in 1:4){
    text(mean(year),sn.data-st,riv[i],font=2)
for(j in 1:n.data[i]){
points(year,ifelse(is.na(fig.dat1[[i]][,j]),NA,sn.data-j-st),col=1,pch=19,cex=1.5)
points(year,ifelse(is.na(fig.dat2[[i]][,j]),NA,sn.data-j-st),col='gray',pch=19,cex=1.5)
axis(4,sn.data-j-st,labels=label.dat[[i]][j],las=1)
}
    st <- st+n.data[i]+1
}

mtext("Year", cex = 1.4,side=1,outer=T)


#'------------------------------------------------------------------------------
#  2.3  Testfishery: Calculate CPUE                     
#'------------------------------------------------------------------------------
# Extract testfish data
testf<-kusko.data[substr(names(kusko.data),1,3)=='rpw']
# combine week 8, 9 and 10 and drop week 9 & 10
testf[,8] <- testf[,8]+testf[,9]+testf[,10]
testf <- testf[,-(9:10)]


# Replace NA to mean proporion for each week 
for (i in 1:dim(testf)[2]) {
  testf[is.na(testf[i]),i] <- colMeans(testf,na.rm=T)[i]
}
testf <- testf[,3:8]

# Extract weekly commercial effort data 
ceff <-kusko.data[substr(names(kusko.data),1,3)=='cew']
# combine week 8, 9  and drop week 9
ceff[,6] <- ceff[,6]+ceff[,7]
ceff <- ceff[,-7]

# Extract weekly commercial catch data
ccat <-kusko.data[substr(names(kusko.data),1,3)=='chw']
# combine week 8, 9  and drop 
ccat[,6] <- ccat[,6]+ccat[,7]
ccat <- ccat[,-7]

# Extract weekly commercial rest data
creg <-kusko.data[substr(names(kusko.data),1,3)=='cfw']
# combine week 8, 9  and drop 
creg[,6] <- pmax(creg[,6],creg[,7])
creg <- creg[,-7]

# Calculate CPUE 
cpue <- ccat/ceff
# Unrestricted Catch preriod is creg = 1
ureg <- ifelse(creg==1,1,0)
# Sum cpue of unrestricted mesh periods 
ur <- rowSums(ureg*cpue,na.rm=TRUE)
# Sum proportion of run periods during the unrestricted 
urp <- rowSums(ureg*testf,na.rm=TRUE)

# Restricted Catch 
ureg <- ifelse(creg==2,1,0)
r <- rowSums(ureg*cpue,na.rm=TRUE)
rp <- rowSums(ureg*testf,na.rm=TRUE)

# Monofilament Catch 
ureg <- ifelse(creg==3|creg==5,1,0)
mono <- rowSums(ureg*cpue,na.rm=TRUE)
monop <- rowSums(ureg*testf,na.rm=TRUE)

c.cpue  <- cbind(ur,r,mono)
c.prop  <- cbind(urp,rp,monop)

#'------------------------------------------------------------------------------
#  2.4  Estimate upper Sub harvest for 1976-1989                              
#'------------------------------------------------------------------------------
# Calculate mean Upriver harvest proportion 
p.l <- with(kusko.data,mean(H.Sub.l/H.Sub,na.rm=TRUE))
kusko.data$H.Sub.l <- with(kusko.data, ifelse(is.na(H.Sub.l),round(H.Sub*p.l),H.Sub.l))

# Calculate total below Bethei Catch  
tlcatch <- rowSums(kusko.data[substr(names(kusko.data),1,2)=='H.'],dims = 1,na.rm=TRUE) - kusko.data$H.Sub
# summarize catches of subsistence fisheries occuring above  Bethel 
tucatch <- with(kusko.data,H.Sub-H.Sub.l)

# Calculate observed minimum escapement
minesc <- rowSums(cbind(w.esc,a.esc), na.rm=TRUE, dims = 1)

# Calculate observed minimum run
minrun <- rowSums(cbind(tlcatch,tucatch,w.esc,a.esc), na.rm=TRUE, dims = 1)

# Calculate the number of observed years
ny <- dim(kusko.data)[1]

fyear <- min(kusko.data$Year)
lyear <- max(kusko.data$Year)

############################################################################
#  2.4  Create ADMB data                                                
############################################################################
dat <- list()
dat$fyear <- fyear
dat$lyear <- lyear 
dat$fage <- 4
dat$lage <- 7
dat$h_low <- tlcatch
dat$h_up <- tucatch
dat$minesc <- minesc
dat$minrun <- minrun
# Radio telemetry data 
dat$mr_tr <- kusko.data$tmr
dat$mrsd_tr <- kusko.data$tmr.sd
# Sonar  data 
dat$sonar <- kusko.data$Sonar
dat$sonar_sd <- kusko.data$Sonar.sd
# MR data
dat$mr <- kusko.data$mr
dat$mrsd <- kusko.data$mr.sd
# MR above Aniak
dat$amr <- kusko.data$amr
dat$amrsd <- kusko.data$amr.sd
#MR Holitna
dat$hmr <- kusko.data$hmr
dat$hmrsd <- kusko.data$hmr.sd
dat$w_esc <- t(w.esc)
dat$a_esc <- t(a.esc)
dat$cpue <- t(c.cpue)
dat$testp <- t(c.prop)
dat$efn_h <- efn_H
dat$age_h <- t(age.class.h)
dat$efn_e <- efn_E
dat$age_e <- t(age.class.e)
dat$SDRec <- 0.5
dat$SDma <- 0.5
# These are assumed CV for weir observation 
# 1:Kwethluk, 2:Tuluksuk, 3: Aniak Salmon, 4:George, 5:Koglukluk, 6:Tatluiksuk, 7:Takotna, 8:Pitka Salmon 
dat$cvw <- c(0.15,0.15,0.15,0.15,0.15,0.15,0.15,0.15)
# These are assumed CV for aerial for likelihood 
# 1:Kwethluk, 2:Tulukuk, 3:Kisaralik, 4:Aniak Salmon, 5:Kipchuk, 6:Anik
# 7:Holokuk, 8:Oskawalik, 9:Holitna, 10:Cheeneetnuk, 11:Gagaryah, 12:Pitka
# 13:Bear, 14: Pitka Salmon
dat$cva <- c(0.50,0.50,0.50,0.50,0.50,0.50,0.50,0.50,0.50,0.50,0.50,0.50,0.50,0.50)
# These are assumed CV for fishery cpue: Large, mix, restrict 
dat$cvf <- c(0.20,0.2,0.2)
# These are assumed cv for harvest 
dat$cvh <- c(0.05,0.05)
# This is a conrol vector for weir likelihood 
# 1:Kwethluk, 2:Tuluksuk, 3: Aniak Salmon, 4:George, 5:Koglukluk, 6:Tatluiksuk, 7:Takotna, 8:Pitka Salmon 
dat$wlike <- c(1,1,1,1,1,1,1,1)
# This is the conrol vector aerial for likelihood 
# 1:Kwethluk, 2:Tulukuk, 3:Kisaralik, 4:Aniak Salmon, 5:Kipchuk, 6:Anik
# 7:Holokuk, 8:Oskawalik, 9:Holitna, 10:Cheeneetnuk, 11:Gagaryah, 12:Pitka
# 13:Bear, 14: Pitka Salmon
dat$alike <- c(1,1,1,1,1,1,1,1,1,1,1,1,1,1)
# Fishery likelihood: Large, mix, restrict  
dat$flike <- c(1,1,1)
# Age compoition likelihood: 
# 1,2 Age comp:harvest, escapement, 3,recruitment deviation, 4,5: maturation deviation 
dat$tflike <- c(0.5,0.1,1,1,1)
# 1, Lower MR, 2. Bethel Sonar, 3. Upper MR, 4. Above Aniak, 5. Holitna
dat$rlike <- c(1,1,1,0,0)
# Year of Fishery selectivity 
dat$sy <- c(1980,1993,2010,2014)

dat <- rapply( dat, f=function(x) ifelse(is.na(x),0,x), how="replace" )

#'------------------------------------------------------------------------------
#  2.3  Data Diaginostic                                                         
#'------------------------------------------------------------------------------

#windows(record=TRUE)
#pairs.linear(cbind(w.esc))
#pairs.linear(cbind(a.esc))
#year <- kusko.data$Year
#par(mfrow=c(5,5))
#na <- dim(a.esc)[2]
#for(i in 1:na){
#plot(year,a.esc[,i],pch=19,yaxt='n')
#}
#nw <- dim(w.esc)[2]
#for(i in 1:nw){
#plot(year,w.esc[,i],pch=19,yaxt='n')
#}

#par(mfrow=c(1,1),mar=c(1,1,1,1),oma=c(1,1,1,1))
#for(i in 1:14){
#par(new=TRUE)
#plot(year,a.esc[,i],type='p',pch=19,col='gray',yaxt='n')

#}
#par(new=TRUE)
#for(i in 1:8){
#plot(year,w.esc[,i],type='p',pch=19,col=1,yaxt='n')
#par(new=TRUE)
#}
#par(new=TRUE)
#plot(year,(kusko.data$Sonar-kusko.data$H.Sub+kusko.data$H.Sub.l), pch=19,cex=1.5,col=3,yaxt='n')
#par(new=TRUE)
#plot(year,(mrs$tmr-kusko.data$H.Com-kusko.data$H.Sub), pch=19,cex=1.2,col=4,yaxt='n')
#par(new=TRUE)
#plot(year,mrs$mr, pch=19,cex=1.2,col=6,yaxt='n')

#'==============================================================================
#   Kuskokwim Chinook Salmon Run reconstruction data construction code
#   Create_data.R 
#   Version 2.0 
#   Written by:   T. Hamaazaki 
#   Date: 07/03/2023
#   Modified: 04/10/2024
#   This program creates an input data for Kuskokowim Run Reconstructiton model
#   The program use kusko.data created  
#   data_file_1:  Harvest, Run, Escapement, fishery data
#   data_file_2:  Harvest and Escapement Age Composition data        
#'==============================================================================
library(reshape2)
#'------------------------------------------------------------------------------
make_RR_Model_data <- function(kusko.data){
#'------------------------------------------------------------------------------
# Calculate the number of observed years
ny <- dim(kusko.data)[1]
fyear <- min(kusko.data$Year)
lyear <- max(kusko.data$Year)
#'------------------------------------------------------------------------------
# Extract testfish data: rpw
#'------------------------------------------------------------------------------
testf<-kusko.data[substr(names(kusko.data),1,3)=='rpw']
# combine week 8, 9 and 10 and drop week 9 & 10
testf[,8] <- testf[,8]+testf[,9]+testf[,10]
testf <- testf[,-(9:10)]
# Replace NA to mean proporion for each week 
for (i in 1:dim(testf)[2]) {
  testf[is.na(testf[i]),i] <- colMeans(testf,na.rm=T)[i]
}
testf <- testf[,3:8]
#'------------------------------------------------------------------------------
# Extract weekly commercial effort and catch data: cew, chw and calculate CPUE 
#'------------------------------------------------------------------------------
ceff <-kusko.data[substr(names(kusko.data),1,3)=='cew']
# combine week 8, 9  and drop week 9
ceff[,6] <- ceff[,6]+ceff[,7]
ceff <- ceff[,-7]

# Extract weekly commercial catch data
ccat <-kusko.data[substr(names(kusko.data),1,3)=='chw']
# combine week 8, 9  and drop 
ccat[,6] <- ccat[,6]+ccat[,7]
ccat <- ccat[,-7]

# Calculate CPUE 
cpue <- ccat/ceff
#'------------------------------------------------------------------------------
# Extract weekly commercial fishery regulation data: cfw 
#'------------------------------------------------------------------------------
creg <-kusko.data[substr(names(kusko.data),1,3)=='cfw']
# combine week 8, 9  and drop 
creg[,6] <- pmax(creg[,6],creg[,7])
creg <- creg[,-7]

# Unrestricted Catch preriod: creg = 1
ureg <- ifelse(creg==1,1,0)
# Sum cpue of unrestricted mesh periods for each year
ur <- rowSums(ureg*cpue,na.rm=TRUE)
# Sum proportion of run periods during the unrestricted 
urp <- rowSums(ureg*testf,na.rm=TRUE)

# Restricted Catch: creg = 2
ureg <- ifelse(creg==2,1,0)
# Sum cpue of restricted mesh periods for each year
r <- rowSums(ureg*cpue,na.rm=TRUE)
# Sum proportion of run periods during the restricted 
rp <- rowSums(ureg*testf,na.rm=TRUE)

# Monofilament Catch: creg = 3 or 5 
ureg <- ifelse(creg==3|creg==5,1,0)
# Sum cpue of monofilament mesh periods for each year
mono <- rowSums(ureg*cpue,na.rm=TRUE)
# Sum proportion of run periods during the monofilament 
monop <- rowSums(ureg*testf,na.rm=TRUE)

c.cpue  <- cbind(ur,r,mono)
c.prop  <- cbind(urp,rp,monop)

#'------------------------------------------------------------------------------
#  2.4  Harvest data                             
#'------------------------------------------------------------------------------
# Calculate mean Upriver harvest proportion 
p.l <- with(kusko.data,mean(H.Sub.l/H.Sub,na.rm=TRUE))
# Apply p.l to missing years
kusko.data$H.Sub.l <- with(kusko.data, ifelse(is.na(H.Sub.l),round(H.Sub*p.l),H.Sub.l))
# Calculate below Bethel Catch  
h_low <- rowSums(kusko.data[substr(names(kusko.data),1,2)=='H.'],dims = 1,na.rm=TRUE) - kusko.data$H.Sub
# Calculate above Bethel Catch
h_up <- with(kusko.data,H.Sub-H.Sub.l)
# Calculate Total catch 
tcatch <-  rowSums(kusko.data[,c('H.Com','H.Sub','H.Sports','H.Test')],dims = 1,na.rm=TRUE)
#'------------------------------------------------------------------------------
#  2.5  Escapement data                             
#'------------------------------------------------------------------------------
# Extract Weier Escapement 
w.esc <- kusko.data[substr(names(kusko.data),1,2)=='w.']
# Extract Aerial Escapement 
a.esc <- kusko.data[substr(names(kusko.data),1,2)=='a.']
# Calculate observed minimum escapement
minesc <- rowSums(cbind(w.esc,a.esc), na.rm=TRUE, dims = 1)
# Calculate observed minimum run
minrun <- rowSums(cbind(tcatch,w.esc,a.esc), na.rm=TRUE, dims = 1)
#'------------------------------------------------------------------------------
#  2.6  Create RR data                             
#'------------------------------------------------------------------------------
dat <- list()
# First Year 
dat$fyear <- fyear
# Last Year 
dat$lyear <- lyear 
# Total Catch 
dat$tcatch <- tcatch
# Catch below Bethel 
dat$h_low <- h_low
# Catth above Bethel
dat$h_up <- h_up
# Mininum Escapement
dat$minesc <- minesc
# Minimum Observed Run
dat$minrun <- minrun
# Inriver Run:  Total run size 
dat$inriv <- kusko.data$In.river
dat$inriv_sd <- kusko.data$In.river.sd
# Lower River Radio telemetry MR  
dat$mr_tr <- kusko.data$tmr
dat$mrsd_tr <- kusko.data$tmr.sd
# Bethel Sonar MR 
dat$sonar <- kusko.data$sonar
dat$sonar_sd <- kusko.data$sonar.sd
# Birch Tree Fishwheel MR 
dat$mr <- kusko.data$mr
dat$mrsd <- kusko.data$mr.sd
# MR above Aniak
dat$amr <- kusko.data$amr
dat$amrsd <- kusko.data$amr.sd
#MR Holitna
dat$hmr <- kusko.data$hmr
dat$hmrsd <- kusko.data$hmr.sd
# Weir Escapement 
dat$w_esc <- t(w.esc)
# Aerial Escapement 
dat$a_esc <- t(a.esc)
dat$cpue <- t(c.cpue)
dat$testp <- t(c.prop)
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
# 1, Lower MR, 2. Bethel Sonar, 3. Upper MR, 4. Above Aniak, 5. Holitna
dat$rlike <- c(1,1,1,0,0)
# Replace NA to 0
dat <- rapply( dat, f=function(x) ifelse(is.na(x),0,x), how="replace" )
return(dat)
}




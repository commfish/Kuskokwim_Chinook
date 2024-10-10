#'==============================================================================
#   Kuskokwim Chinook Salmon Run reconstruction data construction code
#   Create_data.R 
#   Version 2.0 
#   Written by:   T. Hamaazaki 
#   Date: 07/03/2023
#   Modified: 04/10/2024
#   This program creates an input data for Kuskokowim Run Reconstructiton model
#   The program use two data sources 
#   data_file_1:  Harvest, Run, Escapement, fishery data
#   data_file_2:  Harvest and Escapement Age Composition data        
#'==============================================================================

#'------------------------------------------------------------------------------
#   1.0 Source file locations: 
#'------------------------------------------------------------------------------
library(reshape2)
# Set Base directory 
#'------------------------------------------------------------------------------
#  2.1  Read Escapement and Harvest CSV Data                                      
#'------------------------------------------------------------------------------
kusko.data <- read.csv(file.path(Data_dir,data_file.1),header=TRUE, na.string='')

#'------------------------------------------------------------------------------
#  2.2  Testfishery: Calculate CPUE                     
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
#  2.4  Harvest data                             
#'------------------------------------------------------------------------------
# Calculate mean Upriver harvest proportion 
p.l <- with(kusko.data,mean(H.Sub.l/H.Sub,na.rm=TRUE))
# Apply p.l to missing years
kusko.data$H.Sub.l <- with(kusko.data, ifelse(is.na(H.Sub.l),round(H.Sub*p.l),H.Sub.l))
# Calculate total below Bethei Catch  
tlcatch <- rowSums(kusko.data[substr(names(kusko.data),1,2)=='H.'],dims = 1,na.rm=TRUE) - kusko.data$H.Sub
# summarize catches of subsistence fisheries occurring above  Bethel 
tucatch <- with(kusko.data,H.Sub-H.Sub.l)
# calculate Total catch 
t.catch <-  rowSums(kusko.data[,c('H.Com','H.Sub','H.Sports','H.Test')],dims = 1,na.rm=TRUE)

w.esc <- kusko.data[substr(names(kusko.data),1,2)=='w.']
a.esc <- kusko.data[substr(names(kusko.data),1,2)=='a.']


# Calculate observed minimum escapement
minesc <- rowSums(cbind(w.esc,a.esc), na.rm=TRUE, dims = 1)

# Calculate observed minimum run
minrun <- rowSums(cbind(t.catch,w.esc,a.esc), na.rm=TRUE, dims = 1)

# Calculate the number of observed years
ny <- dim(kusko.data)[1]

fyear <- min(kusko.data$Year)
lyear <- max(kusko.data$Year)




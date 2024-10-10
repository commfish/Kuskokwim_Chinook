############################################################################
#  2.4  create RR list data                                                
############################################################################
dat.age <- list()
dat.age$efn_h <- efn_H
dat.age$age_h <- t(age.class.h)
dat.age$efn_e <- efn_E
dat.age$age_e <- t(age.class.e)
dat.age$SDRec <- 0.5  # Assumed Recruit SD
dat.age$SDma <- 0.5   # Assumed Maturity SD
# Age compoition likelihood: 
# 1,2 Age comp:harvest, escapement, 3,recruitment deviation, 4,5: maturation deviation 
dat.age$tflike <- c(0.5,0.1,1,1,1)
dat.age <- rapply(dat.age, f=function(x) ifelse(is.na(x),0,x), how="replace" )


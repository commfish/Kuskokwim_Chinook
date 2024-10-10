#==========================================================================
#  Chinook_RR_Model_default.r
#  This model is the RTMB version of the original ADMB
#  Kuskokwim Run Reconstruction model.  
#  Latest data include Pitka fork weir, but this model excludes the weir 
#  Written by: Toshihide "Hamachan" Hamazaki
#  Date:  04/30/2024
#==========================================================================
Chinook_RR_Model_default <- function(params){
    getAll(TMB.data, params)
#------------------------------------------------------------------------------
#  1.0  Transform parameters  
#------------------------------------------------------------------------------
  t_run <- exp(log_trun)  # Total run
  wesc <- exp(log_wesc)     # slope for weir model
  aesc <- exp(log_aesc)     # slope for aerial model
  q <- exp(log_q)             # slope for catchability
  sd2 <- sqrt(log((exp(log_cvw))^2+1))  #SD parameters Weir
  sd3 <- sqrt(log((exp(log_cva))^2+1))  #SD parameters Aerial
  sd4 <- sqrt(log((exp(log_cvq))^2+1))  #SD parameters cpue
# get imensions 
  nyear <- lyear-fyear+1
  newier <- dim(w_esc)[1]-1  # Pitka fork weir not included. 
  newier <- dim(a_esc)[1]
# Set temprary vector and initialize to 0
  esc <- rep(0,nyear)
  tfw <- rep(0,nweir)  # Weir likelihood
  tfa <- rep(0,naerial) # Aerial likelihood
  tfc <- rep(0,2)  # Com catch likelihood
  tfr <- 0
  f <- 0 

#==============================================================================
# 2.0 Likeihood  
#==============================================================================
 for(i in 1:nyear)
   {
     esc[i] <- t_run[i]-tcatch[i];	   
#===  Total Run ==============================================================     	   
# Total run 
  if(inriv[i]>0)
     {
    sd1 <- sqrt(log((inriv_sd[i]/inriv[i])^2+1))
    tfr <- tfr - dnorm(log(inriv[i]),log(t_run[i]),sd1,TRUE)  
      }   
#============= Escapement =====================================================
# Weir escapement
  for(j in 1:nweir)
    {
	if(w_esc[j,i]>0) 
     {	 
      tfw[j] <- tfw[j] - dnorm(log(w_esc[j,i]),log(esc[i]/wesc[j]),sd2,TRUE) 
     } 
    }
  for(j in 1:naerial)
  {
# Aerial escapement 
	if(a_esc[j,i]>0) 
     {
      tfa[j] <- tfa[j] - dnorm(log(a_esc[j,i]),log(esc[i]/aesc[j]),sd3,TRUE) 
     }  	 
   }     
#'===  CPUE ====================================================================   
	if(cpue[1,i]>0)  
	{
	tfc[1] <- tfc[1] - dnorm(log(cpue[1,i]/testp[1,i]),log(q[1]*t_run[i]),sd4,TRUE)
	}
	if(cpue[2,i]>0)  
	{
	tfc[2] <- tfc[2] - dnorm(log(cpue[2,i]/testp[2,i]),log(q[2]*t_run[i]),sd4,TRUE)
	}
  }           
 
# Sum all likelihood ==========================================================
  f <-  tfr+sum(tfw)+sum(tfa)+sum(tfc)    
#'------------------------------------------------------------------------------
# 3.0   REPORT_SECTION
#'------------------------------------------------------------------------------
  ADREPORT(t_run)
  ADREPORT(esc)
  REPORT(f)
  REPORT(tfw)
  REPORT(tfa)
  REPORT(tfc)
  REPORT(tfr)
  return(f)
  }


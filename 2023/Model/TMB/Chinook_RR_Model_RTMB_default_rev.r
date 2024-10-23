#==========================================================================
#  Chinook_RR_Model_default_rev.r
#  Kuskokwim Run Reconstruction model.   
#  This model keeps original model format and includes Pitka fork weir 
#  Sonar data included
#  Harvest changed to estimates 
#  Written by: Toshihide "Hamachan" Hamazaki
#  Date:  04/30/2024
#==========================================================================
Chinook_RR_Model_RTMB_default_rev <- function(params){
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
  sd5 <- sqrt(log((cvh[1])^2+1))  #SD parameters catch 
#  sd5 <- sqrt(log((exp(log_cvh))^2+1))  #SD parameters catch
# get dimensions 
  nyear <- lyear-fyear+1
  newier <- dim(w_esc)[1] 
  newier <- dim(a_esc)[1]
# Set temprary vector and initialize to 0
  esc <- rep(0,nyear)
  eh_low <- rep(0,nyear)
  eh_up <- rep(0,nyear)  
  low_run <- rep(0,nyear)
  tfw <- rep(0,nweir)  # Weir likelihood
  tfa <- rep(0,naerial) # Aerial likelihood
  tfc <- rep(0,2)  # Com catch likelihood
  tfr <- rep(0,2)
  tfh <- rep(0,2)
  f <- 0 
#==============================================================================
# 2.0 Likeihood  
#==============================================================================
 for(i in 1:nyear)
   {
   eh_low[i] <- exp(-log_fl[i])*t_run[i]
# Run at Bethel is total run minus lower river harvest    
   low_run[i] <- t_run[i] - eh_low[i] 
# Harvest above Bethel   
   eh_up[i] <- exp(-log_fu[i])*low_run[i]   
   esc[i] <- low_run[i] - eh_up[i]
	 
#===  Total Run ==============================================================     	   
# Total run 
  if(inriv[i]>0)
     {
    sd1 <- sqrt(log((inriv_sd[i]/inriv[i])^2+1))
    tfr[1] <- tfr[1] - dnorm(log(inriv[i]),log(t_run[i]),sd1,TRUE)  
      }   
	  
#===  Bethel Sonar  ============================================================     	   
# Sonar
  if(sonar[i]>0)
     { 
    sd6 <- sqrt(log((sonar_sd[i]/sonar[i])^2+1))
    tfr[2] <- tfr[2] - dnorm(log(sonar[i]),log(low_run[i]),sd6,TRUE)  
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
	if(cpue[3,i]>0)  
	{
	tfc[2] <- tfc[2] - dnorm(log(cpue[3,i]/testp[3,i]),log(q[2]*t_run[i]),sd4,TRUE)
	}	
  }           
 
#'==============================================================================   
# Lower River  
  tfh[1] <-  -sum(dnorm(log(h_low),log(eh_low),sd5,TRUE))
# Upriver Harvest
  tfh[2] <- -sum(dnorm(log(h_up),log(eh_up),sd5,TRUE))  

# Sum all likelihood ==========================================================
  f <-  sum(tfr)+sum(tfw)+sum(tfa)+sum(tfc)+sum(tfh)    
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
  REPORT(tfh)
  return(f)
  }


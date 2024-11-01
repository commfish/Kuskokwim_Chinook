#'==============================================================================
#  Chinook_RR_Model_RTMB_New.r
#  This model is the RTMB version of the 
#  Revised Kuksokwim Run Reconstruction model.   
#  Written by: Toshihide "Hamachan" Hamazaki
#  Date: 10/10/2024
#'==============================================================================
Chinook_RR_Model_RTMB_New <- function(params){
    getAll(TMB.data, params)
#'------------------------------------------------------------------------------
#  1.0  Transformed parrameters  
#'------------------------------------------------------------------------------
  r_run <- exp(log_trun)	# Mainstem River run
  wesc <- exp(log_wesc)     # slope for weir model
  aesc <- exp(log_aesc)     # slope for aerial model
  q <- exp(log_q)             # slope for catchability
  sd3 <- sqrt(log((exp(log_cvw))^2+1))  #SD parameters Weir
  sd4 <- sqrt(log((exp(log_cva))^2+1))  #SD parameters Aerial
  sd5 <- sqrt(log((exp(log_cvq))^2+1))  #SD parameters cpue
  sd6 <- sqrt(log((cvh[1])^2+1))  #SD parameters catch

#== Run and Escapment ===========================================  
  esc <- rep(0,nyear)       # Mainstem River escapement 
  low_run <- rep(0,nyear)   # Lower river run
  t_esc <- rep(0,nyear)     # Total escapement Eek included
  t_run <- rep(0,nyear)     # Ttoa Run: Eek included 
#== Harvests ==============================================================
  eh_low <- rep(0,nyear)       # Lower river harvest
  eh_up <- rep(0,nyear)        # Upriver harvest

#==== Likelihood Parameters ===============================================
# Set temprary vector and initialize to 0
  tfw <- rep(0,nweir)  # Weir likelihood
  tfa <- rep(0,naerial) # Aerial likelihood
  tfc <- rep(0,3)  # Com catch likelihood
  tfh <- rep(0,2)	# Harvest likelihood 
  tfr <- rep(0,5)  # Run likelihood
  f <- 0 		   # Likelihood

#------------------------------------------------------------------------------
#  2.0 model_harvest_passage_escapement
#------------------------------------------------------------------------------
 for(i in 1:nyear)
   {
# Lower River harvest 
  eh_low[i] <- exp(-log_fl[i])*r_run[i]
# Run at Bethel is total run minus lower river harvest    
  low_run[i] <- r_run[i] - eh_low[i]
# Harvest above Bethel   
   eh_up[i] <- exp(-log_fu[i])*low_run[i]   
# Drainagewide Escapement is Run at Bethel minus above Behtel harvest 
   esc[i] <-  low_run[i] - eh_up[i] 
# Eek escapement expansion:  Add Eek as Kwethluk Esc x 0.534 
   t_esc[i]  <- esc[i]*(1 + 0.534/wesc[1])
   t_run[i]  <- r_run[i]+ 0.534*esc[i]/wesc[1]   

#'------------------------------------------------------------------------------
#  3.0 Likelidhood 
#'------------------------------------------------------------------------------
# Lower River MR 
  if(mr_tr[i]>0)
     {
    sd1 <- sqrt(log((mrsd_tr[i]/mr_tr[i])^2+1))
    tfr[1] <- tfr[1] - dnorm(log(mr_tr[i]),log(r_run[i]),sd1,TRUE)  
      }   
# Bethel Sonar   	   
   if(sonar[i] >0)
     {
    sd2 <- sqrt(log((sonar_sd[i]/sonar[i])^2+1)); 
	tfr[2] <- tfr[2] - dnorm(log(sonar[i]),log(low_run[i]), sd2,TRUE) 
	 }   
# Kalskage/Birch Tree MR  	   
   if(mr[i] >0)
     {
	sd7 <- sqrt(log((mrsd[i]/mr[i])^2+1))	 
    tfr[3] <- tfr[3]- dnorm(log(mr[i]),log((1-mq)*esc[i]),sd7,TRUE)  
     }   

#============= Above Bethel  Escapement =======================================
# Weir Escapement  
  for(j in 1:nweir)
   {
	if(w_esc[j,i]>0) 
     { 
    tfw[j] <- tfw[j] - dnorm(log(w_esc[j,i]),log(esc[i]/wesc[j]),sd3,TRUE)  
     } 
   }
 
# Aerial Escapements 
  for(l in 1:naerial)
   {
  if(a_esc[l,i]>0) 
     {
    tfa[l] <- tfa[l] - dnorm(log(a_esc[l,i]),log(esc[i]/aesc[l]),sd4,TRUE)  		
     } 
   }   

#'===  CPUE ====================================================================
#' CPUE[2](Restricted period) was removed.   
  if(cpue[1,i]>0)  
	 {
	tfc[1] <- tfc[1] - dnorm(log(cpue[1,i]/testp[1,i]),log(q[1]*r_run[i]),sd5,TRUE)
	 }
  if(cpue[3,i]>0)  
	 {
	tfc[2] <- tfc[2] - dnorm(log(cpue[3,i]/testp[3,i]),log(q[2]*r_run[i]),sd5,TRUE)
	 }	
  }           
 
#'==== Harvest =================================================================   
# Lower River  
  tfh[1] <- -sum(dnorm(log(h_low),log(eh_low),sd6,TRUE))
# Upriver Harvest
  tfh[2] <- -sum(dnorm(log(h_up),log(eh_up),sd6,TRUE))  
   
# Sum all likelihood ===========================================================
  f <- sum(tfw)+sum(tfa)+sum(tfr)+sum(tfc)+sum(tfh);    

#------------------------------------------------------------------------------
# 4.0   REPORT_SECTION
#------------------------------------------------------------------------------
  ADREPORT(t_run);
  ADREPORT(t_esc);
  REPORT(r_run);
  REPORT(esc);
  REPORT(low_run);
  REPORT(eh_up);
  REPORT(f);
  REPORT(tfw);
  REPORT(tfa);
  REPORT(tfc);
  REPORT(tfr);
  REPORT(tfh);
  return (f);
  }


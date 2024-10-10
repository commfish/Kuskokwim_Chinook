#'==============================================================================
#  Chinook_RR_RTMB_New.r
#  This model is the RTMB version of the 
#  Revised Kuksokwim Run Reconstruction model.   
#  Written by: Toshihide "Hamachan" Hamazaki
#  Date:  06/14/2023
#'==============================================================================
Chinook_RR_RTMB_New <- function(params){
    getAll(TMB.data, params)
//------------------------------------------------------------------------------
// 2.0  Define parameters 
//------------------------------------------------------------------------------
  PARAMETER_VECTOR(log_trun);  	//log drainage-wise run
  PARAMETER_VECTOR(log_lwesc);     //log slope for weir Lower
  PARAMETER_VECTOR(log_uwesc);     //log slope for weir Upper
  PARAMETER_VECTOR(log_laesc);	//log slope for aerial Lower
  PARAMETER_VECTOR(log_uaesc); 	//log slope for aerial Upper
  PARAMETER_VECTOR(log_rlwesc);  //log addvar for weir Lower
  PARAMETER_VECTOR(log_ruwesc);  //log addvar for weir Upper
  PARAMETER_VECTOR(log_rlaesc);	//log addvar for aerial Lower
  PARAMETER_VECTOR(log_ruaesc); //log addvar for aerial Upper
  PARAMETER_VECTOR(log_q);      // log catchability  model1
  PARAMETER_VECTOR(log_rq);      // log sd cpue model1
  PARAMETER_VECTOR(log_fl); // annual harvest rate Lower 
  PARAMETER_VECTOR(log_fu); // annual harvest rate Upper
  PARAMETER(mq);          // proporion of lower river esc 
  PARAMETER(sq);          // Sonar fraction 
#'------------------------------------------------------------------------------
#  2.1  Transformed parrameters  
#'------------------------------------------------------------------------------
  t_run <- exp(log_trun)	#Total run
  lwesc <- exp(log_lwesc)     # slope for weir model
  uwesc <- exp(log_uwesc)     # slope for weir model
  laesc <- exp(log_laesc)     # slope for aerial model
  uaesc <- exp(log_uaesc)     # slope for aerial model
  rlwesc <- exp(log_rlwesc)       # slope for weir model
  ruwesc <- exp(log_ruwesc)       # slope for weir model
  rlaesc <- exp(log_rlaesc)       # slope for aerial model
  ruaesc <- exp(log_ruaesc)       # slope for aerial model
  q <- exp(log_q)             # slope for catchability
  rq <- exp(log_rq)           # SD for catchability  
	

#== Run and Escapment ===========================================  
  t_esc <- rep(0,nyear)     # Total escapement 
  low_run <- rep(0,nyear)   # Lower river run
  low_esc <- rep(0,nyear)   # Lower river escapement 
   up_esc <- rep(0,nyear)    # Above Betheil escapement
  
#== Harvests ==============================================================
  eh_low <- rep(0,nyear)       # Lower river harvest
  eh_up <- rep(0,nyear)        # Upriver harvest

#==== Likelihood Parameters ===============================================
# Set temprary vector and initialize to 0
  tfw <- rep(0,8)  # Weir likelihood
  tfa <- rep(0,14) # Aerial likelihood
  tfc <- rep(0,3)  # Com catch likelihood
  tfh <- rep(0,2)	# Harvest likelihood 
  tfr <- rep(0,5)  # Run likelihood

#------------------------------------------------------------------------------
#  4.0 model_harvest_passage_escapement
#------------------------------------------------------------------------------
# Lower River Harvest 
 for(i in 1:nyear)
   {
  eh_low[i] <- exp(-log_fl[i])*t_run[i]
# Run at Bethel is total run minus lower river harvest    
  low_run[i] <- t_run[i] - eh_low[i]
# Harvest above Bethel   
   eh_up[i] <- exp(-log_fu[i])*low_run[i]   
# Drainagewide Escapement is Run at Bethel minus above Behtel harvest 
   t_esc[i] <-  low_run[i] - eh_up[i] 
  }
#==============================================================================
# 5.0 Likeihood  
#==============================================================================
  
 for (int i=0; i<nyear; i++)
   {
# ===  Lower Escapements ======================================================     	   
# Total run 
  if(mr_tr(i)>0)
     {
   Type sd1 <- (log(square(mrsd_tr[i]/mr_tr[i])+1));
    tfr(0) -= dnorm(log(mr_tr(i)),log(t_run(i)),sd1,true);  
      }   
# At Bethel Sonar   	   
   if(sonar(i) >0)
     {
    Type sd2 = sqrt(log(square(sonar_sd(i)/sonar(i))+1)); 
	tfr(1) -= dnorm(log(sonar(i)),log(sq*low_run(i)), sd2,true);  
	 }
#============= Above Bethel  Escapement =======================================
# Kwethluk and Tululkak Weir and Aereial Escapement 
  for(int j=0;j<2;j++)
   {
# Kwethluk and Tululkak Weir Escapement 
	if(wl_esc(j,i)>0) 
     {
	  Type sd3 = sqrt(log(square(cvw(j))+1)+square(rlwesc(j))); 	 
      tfw(j) -= dnorm(log(wl_esc(j,i)),log(mq*t_esc(i)/lwesc(j)),sd3,true);  
     } 
   }
  for(int j=0;j<3;j++)
  {
# Aerial escapement 
	if(al_esc(j,i)>0) 
     {
      Type sd4 = sqrt(log(square(cva(j))+1)+square(rlaesc(j)));
      tfa(j) -= dnorm(log(al_esc(j,i)),log(mq*t_esc(i)/laesc(j)),sd4,true);  
     }  	 
   }   

#============= Above Tulusuk Escapement =======================================   
# Birch Tree MR  	   
   if(mr(i) >0)
     {
	 Type sd5 = sqrt(log(square(mrsd(i)/mr(i))+1));	 
     tfr(2) -= dnorm(log(mr(i)),log((1-mq)*t_esc(i)),sd5,true);  
     }   
# Upriver Weir Escapements 
	 for(int k=0; k<6;k++)
	{
	if(wu_esc(k,i)>0) 
     {
	 Type sd6 = sqrt(log(square(cvw(k+2))+1)+square(ruwesc(k)));
        tfw(k+2) -= dnorm(log(wu_esc(k,i)),log((1-mq)*t_esc(i)/uwesc(k)),sd6,true);  
     } 
    }
# Upriver Aerial Escapements 
  for(int l=0;l<11;l++)
   {
  if(au_esc(l,i)>0) 
     {
      Type sd7 = sqrt(log(square(cva(l+3))+1)+square(ruaesc(l)));
       	tfa(l+3) -= dnorm(log(au_esc(l,i)),log((1-mq)*t_esc(i)/uaesc(l)),sd7,true);  		
     } 
   }   

#===  CPUE ====================================================================   
	if(cpue(0,i)>0)  
	{
	 Type sd8 = sqrt(log(square(cvf(0))+1)+square(rq(0)));	
	tfc(0) -= dnorm(log(cpue(0,i)/testp(0,i)),log(q(0)*t_run(i)),sd8,true);
	}
	if(cpue(1,i)>0)  
	{
	 Type sd9 = sqrt(log(square(cvf(1))+1)+square(rq(1)));		
	tfc(1) -= dnorm(log(cpue(1,i)/testp(1,i)),log(q(1)*t_run(i)),sd9,true);
	}
	if(cpue(2,i)>0)  
	{
	 Type sd10 = sqrt(log(square(cvf(2))+1)+square(rq(2)));		
	tfc(2) -= dnorm(log(cpue(2,i)/testp(2,i)),log(q(2)*t_run(i)),sd10,true);
	}
  }           
 
#====  Harvest ================================================================    
# Lower River  
  tfh(0) -= sum(dnorm(log(h_low),log(eh_low),sqrt(log(square(cvh(0))+1)),true)); 
# Upriver Harvest
  tfh(1) -= sum(dnorm(log(h_up),log(eh_up),sqrt(log(square(cvh(1))+1)),true));  
   
# Sum all likelihood ==========================================================
  Type f = sum(vector<Type>(tfw*wlike))+sum(vector<Type>(tfa*alike))+sum(vector<Type>(tfr*rlike))+sum(vector<Type>(tfc*flike))+sum(tfh);    
#------------------------------------------------------------------------------
# 6.0   REPORT_SECTION
#------------------------------------------------------------------------------
  ADREPORT(t_run);
  ADREPORT(t_esc);
  REPORT(low_run);
#  REPORT(low_esc);
  REPORT(eh_up);
  REPORT(f);
  REPORT(tfw);
  REPORT(tfa);
  REPORT(tfc);
  REPORT(tfr);
  REPORT(tfh);
  return f;
  }


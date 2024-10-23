//==========================================================================
//  Chinook_reconst.cpp
//  This model is the TMB version of the 
//  Kuksokwim Run Reconstruction model based on Escapement.   
//  Written by: Toshihide "Hamachan" Hamazaki
//  Date:  06/14/2023
//==========================================================================
#include <TMB.hpp>
// square function 
template<class Type>
Type square(Type x){
  return pow(x,2);
}

// sqrt function 
template<class Type>
Type sqrt(Type x){
  return pow(x,0.5);
}

template<class Type>
Type objective_function<Type>::operator() ()
{
//------------------------------------------------------------------------------
// 1.0  Data Entry 
//------------------------------------------------------------------------------
  DATA_INTEGER(nyear);  // The number of years 
  DATA_INTEGER(nweir);  // The number of years 
  DATA_INTEGER(naerial);  // The number of years  
  DATA_VECTOR(tcatch);	// Sum of all Catches 
// Read Drainage wide total run size data  
  DATA_VECTOR(inriv);		// Total River MR Estimates
  DATA_VECTOR(inriv_sd);		// Total River MR SD  
// Read Weir data 
  DATA_MATRIX(w_esc);  // Weir Escapement
// Read Aerial data 
  DATA_MATRIX(a_esc);  // Aerial  Escapement 
// Read Weekly Commercial  data  
  DATA_MATRIX(cpue);     // CPUE by fishery 
  DATA_MATRIX(testp);   //  prop of run by fishery  

//------------------------------------------------------------------------------
// 2.0  Define parameters 
//------------------------------------------------------------------------------
  PARAMETER_VECTOR(log_trun);  	//log drainage-wise run
  PARAMETER_VECTOR(log_wesc);    //log slope for weir 
  PARAMETER_VECTOR(log_aesc);	//log slope for aerial 
  PARAMETER(log_cvw);  //log cv for weir
  PARAMETER(log_cva);  //log cv for aerial 
  PARAMETER_VECTOR(log_q);      // log catchability  model1
  PARAMETER(log_cvq);      // log sd cpue model1
 
//------------------------------------------------------------------------------
//  2.1  Transformed parrameters  
//------------------------------------------------------------------------------
  vector<Type> t_run=exp(log_trun);		//Total run
  vector<Type> wesc=exp(log_wesc);     // slope for weir model
  vector<Type> aesc=exp(log_aesc);     // slope for aerial model
  vector<Type> q=exp(log_q);             // slope for catchability
  vector<Type> esc(nyear);
  Type sd2 = sqrt(log(square(exp(log_cvw))+1));
  Type sd3 = sqrt(log(square(exp(log_cva))+1));
  Type sd4 = sqrt(log(square(exp(log_cvq))+1));
  
//==== Likelihood Parameters ===============================================
// Set temprary vector and initialize to 0
  vector<Type>  tfw(nweir);  // Weir likelihood
  vector<Type>  tfa(naerial); // Aerial likelihood
  vector<Type>  tfc(2);  // Com catch likelihood
  Type tfr = 0.0;
  tfw.setZero();
  tfa.setZero();
  tfc.setZero();

//==============================================================================
// 3.0 Likeihood  
//==============================================================================
 for (int i=0; i<nyear; i++)
   {
     esc(i)=t_run(i)-tcatch(i);	   
// ===  Total Run ==============================================================     	   
// Total run 
  if(inriv(i)>0)
     {
   Type sd1 = sqrt(log(square(inriv_sd(i)/inriv(i))+1));
    tfr -= dnorm(log(inriv(i)),log(t_run(i)),sd1,true);  
      }   
//============= Escapement =====================================================
// Weir escapement
  for(int j=0;j<nweir;j++)
   {
	if(w_esc(j,i)>0) 
     {	 
      tfw(j) -= dnorm(log(w_esc(j,i)),log(esc(i)/wesc(j)),sd2,true);  
     } 
   }
  for(int j=0;j<naerial;j++)
  {
// Aerial escapement 
	if(a_esc(j,i)>0) 
     {
      tfa(j) -= dnorm(log(a_esc(j,i)),log(esc(i)/aesc(j)),sd3,true);  
     }  	 
   }     
//===  CPUE ====================================================================   
	if(cpue(0,i)>0)  
	{
	tfc(0) -= dnorm(log(cpue(0,i)/testp(0,i)),log(q(0)*t_run(i)),sd4,true);
	}
	if(cpue(2,i)>0)  
	{
	tfc(1) -= dnorm(log(cpue(2,i)/testp(2,i)),log(q(1)*t_run(i)),sd4,true);
	}
  }           
 
// Sum all likelihood ==========================================================
  Type f = tfr+sum(tfw)+sum(tfa)+sum(tfc);    
//------------------------------------------------------------------------------
// 4.0   REPORT_SECTION
//------------------------------------------------------------------------------
  ADREPORT(t_run);
  ADREPORT(esc);
  REPORT(f);
  REPORT(tfw);
  REPORT(tfa);
  REPORT(tfc);
  REPORT(tfr);
  return f;
  }


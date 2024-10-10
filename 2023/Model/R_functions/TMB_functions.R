#'==============================================================================
# TMB_functions.R
# This code includes TMB related functions
#'==============================================================================

#'------------------------------------------------------------------------------
#' runRTMB
#' Reads data, model, and parameter, and run RTMB model 
#'------------------------------------------------------------------------------
runRTMB<-function(rtmb_data,rtmb_fcn,params,upper,lower){
  environment(rtmb_fcn) <- list2env(rtmb_data) 
  #--some function that return a "data" list
  objfun <- RTMB::MakeADFun(rtmb_fcn,params) #--rtmb_fcn is the source'd objective function that uses "data"  
  fit <- nlminb(objfun$par, objfun$fn, control = list(eval.max = 10000, iter.max=10000,trace=1),upper = upper, loewr = lower)  
  return(objfun)#--returns the "made" RTMB function
}
## Show Color Shades -------------------------------------------------------   
tcol <- function(color, percent = 50, name = NULL) {
  #	  color = color name
  #	percent = % transparency
  #	   name = an optional name for the color
  ## Get RGB values for named color
  rgb.val <- col2rgb(color)
  ## Make new color using input color as base and alpha set by transparency
  t.col <- rgb(rgb.val[1], rgb.val[2], rgb.val[3],
               max = 255,
               alpha = (100-percent)*255/100,
               names = name)
  return(t.col)
}

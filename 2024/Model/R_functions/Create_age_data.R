#'------------------------------------------------------------------------------
#  Create_age_data.R
#  This creates age data 
#'------------------------------------------------------------------------------
make_age_data <- function(age.data){
# Extract age data 
age.class <- age.data[!names(age.data) %in% c('n_H','n_E')]
# Change wide to long: This will put column name as "variable" and age comp as "value"
age.class <- melt(age.class, id.vars='Year')
# Extract Havest(H) vs Escapement (E) from the variable column 
age.class$loc <- substr(age.class$variable,1,1)
# Extract age from the variable column
age.class$at <- as.numeric(substr(age.class$variable,3,5))
# Convert European age notation to simple age 
age.class$Age <- with(age.class, floor(at)+10*(at-floor(at))+1)
# Sum proportion by age 
age.class.sum <- aggregate(value~loc+Age+Year, sum, data=age.class)
# Change back long to Wide format 
age.class  <- dcast(age.class.sum,Year~loc+Age,value.var='value')
# Get Harvest age (select column starting from H)
age.class.h <- age.class[substr(names(age.class),1,1)=='H']
# Get Standardized proprtion 
age.class.h <- age.class.h/rowSums(age.class.h)  
# Extract Escapement age (select column starting from E) 
age.class.e <- age.class[substr(names(age.class),1,1)=='E']
# Get Standardized proprion 
age.class.e <- age.class.e/rowSums(age.class.e)  
# Select age sample size for escapement and harvest 
dat <- data.frame(Year=age.data$Year, age.class.h, age.class.e, Efh = age.class$Efn)
return(dat)
}

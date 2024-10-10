##########################################################################################
#   ADMB DATA Write and Report Read Program
#   function makedata creates ADMB readable data file. 
#   dat: R list data  outfn (output file name)
#   function reptRlist  reads ADMB output. 
#   fn: ADMB .ret output file name
##########################################################################################

##########################################################################################
#   Read ADMB dat file to R
##########################################################################################
datatoR = function(fn,par=F)
{
if (par) flush=F else flush=T
ifile=scan(fn,what="character",flush=flush, blank.lines.skip=F,quiet=TRUE)
idx=sapply(as.double(ifile),is.na)
vnam=ifile[idx] #list names
nv=length(vnam) #number of objects
A=list()
ir=0
for(i in 1:nv)
{
	ir=match(vnam[i],ifile)
	if(i!=nv) irr=match(vnam[i+1],ifile) else irr=length(ifile)+1 #next row
	dum=NA
	if(par) dum=as.numeric(as.character(ifile[(ir+1):(irr-1)])) else {
		if(irr-ir==2) dum=as.double(scan(fn,skip=ir,nlines=1,quiet=T,what="")) 
		if(irr-ir>2) dum=as.matrix(read.table(fn,skip=ir,nrow=irr-ir-1,fill=T))
		}

	if(is.numeric(dum))#Logical test to ensure dealing with numbers
	name=sub('#','',vnam[i])
	A[[name]]=dum
}
return(A)
}

makedata <- function(dat,outfn){
# Save original file as filename_0.dat
file.remove(paste(outfn))
out_file <- file(paste(outfn), open='a')
for (j in seq_along(dat)){
    write.table(paste0("#",names(dat)[j]), file=out_file, sep="", dec=".", quote=FALSE, col.names=FALSE, row.names=FALSE)  
	write(t(dat[[j]]), file=out_file, sep=" ", 
	ncolumns= ifelse(is.numeric(dim(dat[[j]])[2]),dim(dat[[j]])[2],length(dat[[j]])))
	}
close(out_file) 
}

# Define the report data
reptoRlist = function(fn,par=F)
{
if (par) flush=F else flush=T
ifile=scan(fn,what="character",flush=flush, blank.lines.skip=F,quiet=T)
idx=sapply(as.double(ifile),is.na)
vnam=ifile[idx] #list names
nv=length(vnam) #number of objects
A=list()
ir=0
for(i in 1:nv)
{
	ir=match(vnam[i],ifile)
	if(i!=nv) irr=match(vnam[i+1],ifile) else irr=length(ifile)+1 #next row
	dum=NA
	if(par) dum=as.numeric(as.character(ifile[(ir+1):(irr-1)])) else {
		if(irr-ir==2) dum=as.double(scan(fn,skip=ir,nlines=1,quiet=T,what="")) 
		if(irr-ir>2) dum=as.matrix(read.table(fn,skip=ir,nrow=irr-ir-1,fill=T))
		}

	if(is.numeric(dum))#Logical test to ensure dealing with numbers
	name=sub(':','',vnam[i])
	A[[name]]=dum
}
return(A)
}

getADMBHessian <- function(){
## This function reads in all of the information contained in the
## admodel.hes file. Some of this is needed for relaxing the covariance
## matrix, and others just need to be recorded and rewritten to file so ADMB
## "sees" what it’s expecting.
filename <- file("admodel.hes", "rb")
on.exit(close(filename))
num.pars <- readBin(filename, "integer", 1)
hes.vec <- readBin(filename, "numeric", num.pars^2)
hes <- matrix(hes.vec, ncol=num.pars, nrow=num.pars)
hybrid_bounded_flag <- readBin(filename, "integer", 1)
scale <- readBin(filename, "numeric", num.pars)
result <- list(num.pars=num.pars, hes=hes,
hybrid_bounded_flag=hybrid_bounded_flag, scale=scale)
return(result)
}


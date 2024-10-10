#'==============================================================================
#  Create data summary figures
#  This code reads lookup data and create available data summary plot
#'==============================================================================
# Read lookup table
lookup <- read.csv(file.path(Data_dir,data_file.3),header = TRUE, na.string='')
# Extract Weir escapement data
# Mark-recapcutre data 
mr <- lookup[which(lookup$River =='MR'),]
mrs <- kusko.data[,names(kusko.data) %in% mr$CSV]


# Extract Upriver data 
up <- lookup[which(lookup$River =='Up'),]
Esc.fig.u1 <- kusko.data[,names(kusko.data) %in% up$CSV]
Esc.fig.u2 <- kusko.data[,names(kusko.data) %in% up$CSV]
mid <- lookup[which(lookup$River =='Mid'),]
Esc.fig.m1 <- kusko.data[,names(kusko.data) %in% mid$CSV]
Esc.fig.m2 <- kusko.data[,names(kusko.data) %in% mid$CSV]
low <- lookup[which(lookup$River =='Low'),]
Esc.fig.l1 <- kusko.data[,names(kusko.data) %in% low$CSV]
Esc.fig.l2 <- kusko.data[,names(kusko.data) %in% low$CSV]
mr <- lookup[which(lookup$River =='MR'),]
M.fig.s1 <- cbind(kusko.data[,names(kusko.data) %in% mr$CSV],kusko.data$sonar)
M.fig.s2 <- cbind(mrs,kusko.data$sonar)
Cpue.fig.s1 <- c.cpue
Cpue.fig.s1[Cpue.fig.s1==0]<- NA
Cpue.fig.s2 <- c.cpue
# Erase data 
Esc.fig.u2[,] <- NA
Esc.fig.m2[,] <- NA
Esc.fig.l2[,] <- NA
M.fig.s2[,] <-NA
Cpue.fig.s2[,] <- NA
year <- kusko.data$Year
fig.dat1 <- list(Esc.fig.u1,Esc.fig.m1,Esc.fig.l1,M.fig.s1,Cpue.fig.s1)
fig.dat2 <- list(Esc.fig.u2,Esc.fig.m2,Esc.fig.l2,M.fig.s2,Cpue.fig.s2)
n.u <- dim(Esc.fig.u1)[2]
n.m <- dim(Esc.fig.m1)[2]
n.l <- dim(Esc.fig.l1)[2]
n.s <- dim(M.fig.s1)[2]
n.c <- dim(Cpue.fig.s1)[2]

n.data <- c(n.u, n.m, n.l, n.s,n.c)
sn.data <- sum(n.data)+4
riv <- c('Upper','Middle','Lower','Mainstem','CPUE')
label.u <- up$Label
label.m <- mid$Label
label.l <- low$Label	 
label.s <- c(mr$Label,'Bethel Sonar')	
label.c <- c('Unrestricted','Restricted','Monofilament')	
label.dat <- list(label.u,label.m,label.l,label.s,label.c)  

par(mfrow=c(1,1),mar=c(2,2,2,9),oma=c(3,1,1,1))
plot(year,Esc.fig.u1[,1],xlim=c(min(year)-1,max(year+1)),ylim=c(0,sn.data),xaxs='i',yaxt='n',
       type="n",xlab="",ylab="",main="",cex.main=1.)  
st <- 0
for (i in 1:5){
    text(mean(year),sn.data-st,riv[i],font=2)
for(j in 1:n.data[i]){
points(year,ifelse(is.na(fig.dat1[[i]][,j]),NA,sn.data-j-st),col=1,pch=19,cex=1.5)
points(year,ifelse(is.na(fig.dat2[[i]][,j]),NA,sn.data-j-st),col='gray',pch=19,cex=1.5)
axis(4,sn.data-j-st,labels=label.dat[[i]][j],las=1)
}
    st <- st+n.data[i]+1
}

mtext("Year", cex = 1.4,side=1,outer=T)


library(readxl)
Survival_example_mean <- read_excel("C:/Users/dell/Desktop/2023.0257/scr/output/Survival_example10_mean.xlsx",col_names = FALSE) #mean
#Survival_example_std <- read_excel("C:/Users/dell/Desktop/2023.0257/scr/output/Survival_example10_std.xlsx",col_names = FALSE) # standard deviation


nme=13                                           
xx=c('50','100','200');                          
ns=c('T-PL','A-PL')
nn=length(ns)
awei=array(0,dim=c(nn,length(xx),nme))

# mean
awei[1,,]=as.matrix(Survival_example_mean[1:3,1:13],3,13) 
awei[2,,]=as.matrix(Survival_example_mean[4:6,1:13],3,13) 

# standard deviation
#awei[1,,]=as.matrix(Survival_example_std[1:3,1:13],3,13) 
#awei[2,,]=as.matrix(Survival_example_std[4:6,1:13],3,13) 


lthick1=matrix(2,nme,1)
col1=c("green","red","orange","pink","gray","steelblue","lightblue","maroon","purple","blue","black","brown","darkgreen")
ltype1=c(3,10,9,8,7,4,5,8,1,6,1,1,2)
par(mfrow=c(1,2))
for (ni in 1:nn) {
  wei=awei[ni,,]
  par(mar=c(5, 5, 2, 1))
  if (ni==1) {
    # mean
    matplot(wei,type="l",xlab=expression(n),ylab="Relative Targeted Empirical Ranking Risk",lty=ltype1,lwd=lthick1,ylim=c(0.95,1.2),xaxt='n',col=col1,mgp=c(2.2,1,0))
    matpoints(wei,type="p",pch=c(0,1,9:14,2,3,4,5,6),main=main1,xlab=expression(c),ylab="Relative Targeted Empirical Ranking Risk",lty=ltype1,lwd=lthick1,ylim=c(0.8,2),xaxt='n',col=col1,mgp=c(2.2,1,0))
    # standard deviation
#    matplot(wei,type="l",xlab=expression(n),ylab="Relative Standard Deviation of Targeted Empirical Ranking Risk",lty=ltype1,lwd=lthick1,ylim=c(0.6,1.8),xaxt='n',col=col1,mgp=c(2.2,1,0))
#    matpoints(wei,type="p",pch=c(0,1,9:14,2,3,4,5,6),main=main1,xlab=expression(c),ylab="Relative Targeted Empirical Ranking Risk",lty=ltype1,lwd=lthick1,ylim=c(0.8,2),xaxt='n',col=col1,mgp=c(2.2,1,0))
  }else{
    # mean
    matplot(wei,type="l",xlab=expression(n),ylab="Relative Surrogate Empirical Ranking Risk",lty=ltype1,lwd=lthick1,ylim=c(0.98,1.03),xaxt='n',col=col1,mgp=c(2.2,1,0))
    matpoints(wei,type="p",pch=c(0,1,9:14,2,3,4,5,6),main=main1,xlab=expression(c),ylab="Relative Surrogate Empirical Ranking Risk",lty=ltype1,lwd=lthick1,ylim=c(0.7,3),xaxt='n',col=col1,mgp=c(2.2,1,0))
    # standard deviation
#    matplot(wei,type="l",xlab=expression(n),ylab="Relative Standard Deviation of Surrogate Empirical Ranking Risk",lty=ltype1,lwd=lthick1,ylim=c(0.5,4),xaxt='n',col=col1,mgp=c(2.2,1,0))
#    matpoints(wei,type="p",pch=c(0,1,9:14,2,3,4,5,6),main=main1,xlab=expression(c),ylab="Relative Surrogate Empirical Ranking Risk",lty=ltype1,lwd=lthick1,ylim=c(0.7,3),xaxt='n',col=col1,mgp=c(2.2,1,0))
  }
  axis(1, 1:3, c('50','100','200'))
  legend('topleft',c("Full","AIC","BIC","SAIC","SBIC","MMA","Equal","K-data","RMA-op","RMA-2","RMA-5","RMA-10","RMA-n/2"), lty=ltype1,pch=c(0,1,9:14,2,3,4,5,6),col=col1,bg='white',cex=0.9)
}

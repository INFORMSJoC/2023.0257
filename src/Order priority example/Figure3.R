library(readxl)
OP_example_mean <- read_excel("C:/Users/dell/Desktop/2023.0257/scr/output/OP_example_mean.xlsx",col_names = FALSE) #mean
#OP_example_std <- read_excel("C:/Users/dell/Desktop/2023.0257/scr/output/OP_example_std.xlsx",col_names = FALSE)  #standard deviation

nme=12                                       
xx=c('15','20','25','28');                         
ns=c('T-PL','A-PL')
nn=length(ns)
awei=array(0,dim=c(nn,length(xx),nme))

#mean

awei[1,,]=as.matrix(OP_example_mean[1:4,1:12],4,12) 
awei[2,,]=as.matrix(OP_example_mean[5:8,1:12],4,12) 

#srandard deviation

#awei[1,,]=as.matrix(OP_example_std[1:4,1:12],4,12) 
#awei[2,,]=as.matrix(OP_example_std[5:8,1:12],4,12) 

lthick1=matrix(2,nme,1)
col1=c("green","red","orange","pink","gray","steelblue","lightblue","maroon","purple","blue","black","darkgreen")
ltype1=c(3,10,9,8,7,4,5,8,1,6,1,2)
par(mfrow=c(1,2))
for (ni in 1:nn) {
  wei=awei[ni,,]
  par(mar=c(5, 5, 2, 1))
  if (ni==1) {
    #mean
    matplot(wei,type="l",xlab=expression(n),ylab="Relative Targeted Empirical Ranking Risk",lty=ltype1,lwd=lthick1,ylim=c(0.8,2),xaxt='n',col=col1,mgp=c(2.2,1,0))
    matpoints(wei,type="p",pch=c(0,1,9:14,2,3,4,6),main=main1,xlab=expression(c),ylab="Relative Targeted Empirical Ranking Risk",lty=ltype1,lwd=lthick1,ylim=c(0.8,2.5),xaxt='n',col=col1,mgp=c(2.2,1,0))
    #srandard deviation
#    matplot(wei,type="l",xlab=expression(n),ylab="Relative Standard Deviation of Targeted Empirical Ranking Risk",lty=ltype1,lwd=lthick1,ylim=c(0.96,3.5),xaxt='n',col=col1,mgp=c(2.2,1,0))
#    matpoints(wei,type="p",pch=c(0,1,9:14,2,3,4,6),main=main1,xlab=expression(c),ylab="Relative Standard Deviation of Targeted Empirical Ranking Risk",lty=ltype1,lwd=lthick1,ylim=c(0.8,2.5),xaxt='n',col=col1,mgp=c(2.2,1,0))
    
      }else{
    #mean    
    matplot(wei,type="l",xlab=expression(n),ylab="Relative Surrogate Empirical Ranking Risk",lty=ltype1,lwd=lthick1,ylim=c(0.96,1.34),xaxt='n',col=col1,mgp=c(2.2,1,0))
    matpoints(wei,type="p",pch=c(0,1,9:14,2,3,4,6),main=main1,xlab=expression(c),ylab="Relative Surrogate Empirical Ranking Risk",lty=ltype1,lwd=lthick1,ylim=c(1,1.2),xaxt='n',col=col1,mgp=c(2.2,1,0))
    #srandard
#    matplot(wei,type="l",xlab=expression(n),ylab="Relative Standard Deviation of Surrogate Empirical Ranking Risk",lty=ltype1,lwd=lthick1,ylim=c(0.7,2.5),xaxt='n',col=col1,mgp=c(2.2,1,0))
#    matpoints(wei,type="p",pch=c(0,1,9:14,2,3,4,6),main=main1,xlab=expression(c),ylab="Relative Standard Deviation of Surrogate Empirical Ranking Risk",lty=ltype1,lwd=lthick1,ylim=c(1,1.2),xaxt='n',col=col1,mgp=c(2.2,1,0))
      }
  axis(1, 1:4, c('15','20','25','28'))
  legend('topleft',c("Full","AIC","BIC","SAIC","SBIC","MMA","Equal","K-data","RMA-op","RMA-2","RMA-5","RMA-n/2"), lty=ltype1, pch=c(0,1,9:14,2,3,4,6),col=col1,bg='white',cex=0.9)
}

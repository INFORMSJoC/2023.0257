library(readxl)
Piston_example_mean <- read_excel("C:/Users/dell/Desktop/2023.0257/scr/output/SimA6_Piston_example_mean.xlsx",col_names = FALSE)
#Piston_example_std <- read_excel("C:/Users/dell/Desktop/2023.0257/scr/output/SimA6_Piston_example_std.xlsx",col_names = FALSE)

nme=13                                       
xx=c('15','20','25','28');                         
ns=c('T-PL','A-PL')
nn=length(ns)
awei=array(0,dim=c(nn,length(xx),nme))

# mean
awei[1,,]=as.matrix(Piston_example_mean[1:4,1:13],4,13) 
awei[2,,]=as.matrix(Piston_example_mean[5:8,1:13],4,13) 

# std
#awei[1,,]=as.matrix(Piston_example_std[1:4,1:13],4,13) 
#awei[2,,]=as.matrix(Piston_example_std[5:8,1:13],4,13) 

#sam_size=1000;
#n=c(0.4,0.45,0.5,0.55,0.6,0.65,0.7,0.75,0.8,0.85,0.9,0.95,1)

lthick1=matrix(2,nme,1)
#col1=c("purple","pink","gray","green","red","blue","lightgreen","orange","black","brown","deepred","lightblue")
col1=c("green","red","orange","pink","gray","steelblue","lightblue","maroon","purple","blue","black","brown","darkgreen")
ltype1=c(3,10,9,8,7,4,5,8,1,6,1,1,2)
#layout(matrix(c(1,2,3,4,5), ncol=2, byrow=TRUE), heights=c(4,4,0.5))
par(mfrow=c(1,2))
#pdf(file=figname,width = 6, height = 4)
for (ni in 1:nn) {
  wei=awei[ni,,]
  #  main1 = bquote(paste('n','=',.(ns[ni])))
  par(mar=c(5, 5, 2, 1))
  if (ni==1) {
    matplot(wei,type="l",xlab=expression(n),ylab="Relative Targeted Empirical Ranking Risk",lty=ltype1,lwd=lthick1,ylim=c(0.95,1.2),xaxt='n',col=col1,mgp=c(2.2,1,0))
    matpoints(wei,type="p",pch=c(0,1,9:14,2,3,4,5,6),main=main1,xlab=expression(c),ylab="Relative Targeted Empirical Ranking Risk",lty=ltype1,lwd=lthick1,ylim=c(0.8,2),xaxt='n',col=col1,mgp=c(2.2,1,0))
  }else{
    matplot(wei,type="l",xlab=expression(n),ylab="Relative Surrogate Empirical Ranking Risk",lty=ltype1,lwd=lthick1,ylim=c(0.95,1.3),xaxt='n',col=col1,mgp=c(2.2,1,0))
    matpoints(wei,type="p",pch=c(0,1,9:14,2,3,4,5,6),main=main1,xlab=expression(c),ylab="Relative Surrogate Empirical Ranking Risk",lty=ltype1,lwd=lthick1,ylim=c(0.7,3),xaxt='n',col=col1,mgp=c(2.2,1,0))
  }
  axis(1, 1:3, c('50','100','200'))
  #if(ni==1) 
  legend('topleft',c("Full","AIC","BIC","SAIC","SBIC","MMA","Equal","K-data","RMA-op","RMA-2","RMA-5","RMA-10","RMA-n/2"), lty=ltype1,pch=c(0,1,9:14,2,3,4,5,6),col=col1,bg='white',cex=0.7)
}

#figname <- paste("figure_sample_",sam_size,".eps",sep="")
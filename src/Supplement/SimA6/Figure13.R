library(readxl)
Piston_example_compare_mean <- read_excel("C:/Users/dell/Desktop/2023.0257/scr/output/Piston_example_compare_mean.xlsx",col_names = FALSE)
#Piston_example_compare_std <- read_excel("C:/Users/dell/Desktop/2023.0257/scr/output/Piston_example_compare_std.xlsx",col_names = FALSE)

nme=2                                        
xx=c('40','60','80','90');                         
ns=c('T-PL','A-PL')
nn=length(ns)
awei=array(0,dim=c(nn,length(xx),nme))

# mean
awei[1,,]=as.matrix(Piston_example_compare_mean[1:4,c(1,2)],4,2) 
awei[2,,]=as.matrix(Piston_example_compare_mean[5:8,c(1,2)],4,2) 

# std
#awei[1,,]=as.matrix(Piston_example_compare_std[1:4,c(1,2)],4,2) 
#awei[2,,]=as.matrix(Piston_example_compare_std[5:8,c(1,2)],4,2) 


lthick1=matrix(2,nme,1)
#col1=c("purple","pink","gray","green","red","blue","lightgreen","orange","black","brown","deepred","lightblue")
col1=c("blue","red")
ltype1=c(1,1)
#layout(matrix(c(1,2,3,4,5), ncol=2, byrow=TRUE), heights=c(4,4,0.5))
par(mfrow=c(1,2))
#pdf(file=figname,width = 6, height = 4)
for (ni in 1:nn) {
  wei=awei[ni,,]
  #  main1 = bquote(paste('n','=',.(ns[ni]))) Standard Deviation of 
  par(mar=c(5, 5, 2, 1))
  if (ni==1) {
    matplot(wei,type="l",xlab=expression(n),ylab="Relative Targeted Empirical Ranking Risk",lty=ltype1,lwd=lthick1,ylim=c(0.75,1.6),xaxt='n',col=col1,mgp=c(2.2,1,0))
    matpoints(wei,type="p",pch=c(0,1),main=main1,xlab=expression(c),ylab="Relative Targeted Empirical Ranking Risk",lty=ltype1,lwd=lthick1,ylim=c(0.6,1),xaxt='n',col=col1,mgp=c(2.2,1,0))
  }else{
    matplot(wei,type="l",xlab=expression(n),ylab="Relative Surrogate Empirical Ranking Risk",lty=ltype1,lwd=lthick1,ylim=c(0.9,1.2),xaxt='n',col=col1,mgp=c(2.2,1,0))
    matpoints(wei,type="p",pch=c(0,1),main=main1,xlab=expression(c),ylab="Relative Surrogate Empirical Ranking Risk",lty=ltype1,lwd=lthick1,ylim=c(0.7,1),xaxt='n',col=col1,mgp=c(2.2,1,0))
  }
  axis(1, 1:4, c('40','60','80','90'))
  #if(ni==1) 
  legend('topright',c("Correlation","BIC"), lty=ltype1,pch=c(0,1),col=col1,bg='white',cex=0.9)
}
#figname <- paste("figure_sample_",sam_size,".eps",sep="")
#savePlot(file=figname,type="eps",dev.cur())
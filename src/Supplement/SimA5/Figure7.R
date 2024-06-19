library(readxl)
SimA5_compare <- read_excel("C:/Users/dell/Desktop/2023.0257/scr/output/SimA5_compare.xlsx",col_names = FALSE)

nme=3                                                                    
c=c(0.01,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,0.99);                          
ns=c(100,200,500,1000)
nn=length(ns)
awei=array(0,dim=c(nn,length(c),nme))

#Surrogate PL

#awei[1,,]=as.matrix(SimA5_compare[1:11,c(1:3)],11,3) 
#awei[2,,]=as.matrix(SimA5_compare[25:35,c(1:3)],11,3) 
#awei[3,,]=as.matrix(SimA5_compare[49:59,c(1:3)],11,3) 
#awei[4,,]=as.matrix(SimA5_compare[73:83,c(1:3)],11,3)  

#Targeted PL

awei[1,,]=as.matrix(SimA5_compare[12:22,c(1:3)],11,3) 
awei[2,,]=as.matrix(SimA5_compare[36:46,c(1:3)],11,3) 
awei[3,,]=as.matrix(SimA5_compare[60:70,c(1:3)],11,3) 
awei[4,,]=as.matrix(SimA5_compare[84:94,c(1:3)],11,3)  

lthick1=matrix(2,nme,1)
col1=c("green","red","black")
ltype1=c(1,1,1)
ylim1=c(0.9,1.3)
par(mfrow=c(2,2))
#pdf(file=figname,width = 6, height = 4 )
for (ni in 1:nn) {
  wei=awei[ni,,]
  main1 = bquote(paste('n','=',.(ns[ni])))
  par(mar=c(5, 5, 2, 1))
  matplot(c,wei,type="l",main=main1,xlab=expression(c),ylab="Relative Targeted Empirical Ranking Risk",lty=ltype1,lwd=lthick1,ylim=ylim1,xaxt='n',col=col1,mgp=c(2.2,1,0))
  matpoints(c,wei,type="p",pch=0:2,main=main1,xlab=expression(c),ylab="Relative Targeted Empirical Ranking Risk",lty=ltype1,lwd=lthick1,ylim=ylim1,xaxt='n',col=col1,mgp=c(2.2,1,0))
  axis(1,seq(0,1,0.2))
  #if(ni==1) 
  legend('topright',c("RMA-Normal","RMA-Log","RMA-Mix"), lty=ltype1,pch=0:6,col=col1,bg='white',cex=0.7)
}
#figname <- paste("figure_sample_",sam_size,".eps",sep="")
savePlot(file=figname,type="eps",dev.cur())

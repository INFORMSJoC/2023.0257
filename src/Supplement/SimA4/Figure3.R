library(readxl)
SimA4_PL <- read_excel("C:/Users/dell/Desktop/2023.0257/scr/output/SimA4_PL.xlsx",col_names = FALSE)

nme=9                                                                  
cx=c(0.01,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,0.99);                      
ns=c(100,200,500,1000)
nn=length(ns)
awei=array(0,dim=c(nn,length(cx),nme))

#Surrogate PL

#awei[1,,]=as.matrix(SimA4_PL[1:11,c(1:9)],11,9) 
#awei[2,,]=as.matrix(SimA4_PL[25:35,c(1:9)],11,9) 
#awei[3,,]=as.matrix(SimA4_PL[49:59,c(1:9)],11,9) 
#awei[4,,]=as.matrix(SimA4_PL[73:83,c(1:9)],11,9)  

#Targeted PL

awei[1,,]=as.matrix(SimA4_PL[12:22,c(1:9)],11,9) 
awei[2,,]=as.matrix(SimA4_PL[36:46,c(1:9)],11,9) 
awei[3,,]=as.matrix(SimA4_PL[60:70,c(1:9)],11,9) 
awei[4,,]=as.matrix(SimA4_PL[84:94,c(1:9)],11,9)  

lthick1=matrix(2,nme,1)
col1=c("steelblue","lightblue","green","red","purple","blue","black","brown","darkgreen")
ltype1=c(4,5,3,10,1,6,1,1,2)
ylim1=c(0.89,2.5)
par(mfrow=c(2,2))
for (ni in 1:nn) {
  wei=awei[ni,,]
  main1 = bquote(paste('n','=',.(ns[ni])))
  par(mar=c(5, 5, 2, 1))
  #Surrogate PL
#  matplot(cx,wei, type="l",main=main1,xlab=expression(c),ylab="Relative Surrogate Empirical Ranking Risk",lty=ltype1,lwd=lthick1,ylim=ylim1,xaxt='n',col=col1,mgp=c(2.2,1,0))
#  matpoints(cx,wei,type="p",pch=c(12,13,0,1,2,3,4,5,6),main=main1,xlab=expression(c),ylab="Relative Surrogate Empirical Ranking Risk",lty=ltype1,lwd=lthick1,ylim=ylim1,xaxt='n',col=col1,mgp=c(2.2,1,0))
  #Targeted PL
  matplot(cx,wei, type="l",main=main1,xlab=expression(c),ylab="Relative Targeted Empirical Ranking Risk",lty=ltype1,lwd=lthick1,ylim=ylim1,xaxt='n',col=col1,mgp=c(2.2,1,0))
  matpoints(cx,wei,type="p",pch=c(12,13,0,1,2,3,4,5,6),main=main1,xlab=expression(c),ylab="Relative Targeted Empirical Ranking Risk",lty=ltype1,lwd=lthick1,ylim=ylim1,xaxt='n',col=col1,mgp=c(2.2,1,0))
  axis(1,seq(0,1,0.2))
  legend('topleft',c("Incorrect 1","Incorrect 2","Incorrect 3","True","RMA-op","RMA-2","RMA-5","RMA-10","RMA-n/2"), lty=ltype1,pch=c(12,13,0,1,2,3,4,5,6),col=col1,bg='white',cex=0.71)
}

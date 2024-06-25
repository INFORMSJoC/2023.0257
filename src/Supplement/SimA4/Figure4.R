library(readxl)
SimA4_w <- read_excel("C:/Users/dell/Desktop/2023.0257/scr/output/SimA4_w.xlsx",col_names = FALSE)

nme=5                                                                     
c=c(0.01,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,0.99);                        
ns=c(100,200,500,1000)
nn=length(ns)
awei=array(0,dim=c(nn,length(c),nme))

#weight
awei[1,,]=as.matrix(SimA4_w[1:11,c(1:5)],11,5) 
awei[2,,]=as.matrix(SimA4_w[14:24,c(1:5)],11,5) 
awei[3,,]=as.matrix(SimA4_w[27:37,c(1:5)],11,5) 
awei[4,,]=as.matrix(SimA4_w[40:50,c(1:5)],11,5) 


lthick1=matrix(2,nme,1)
col1=c("purple","blue","black","brown","darkgreen")
ltype1=c(1,6,1,1,2)
ylim1=c(0.4,1.01)
par(mfrow=c(2,2))
for (ni in 1:nn) {
  wei=awei[ni,,]
  main1 = bquote(paste('n','=',.(ns[ni])))
  par(mar=c(5, 5, 2, 1))
  matplot(c,wei,type="l",main=main1,xlab=expression(c),ylab="Weight",lty=ltype1,lwd=lthick1,ylim=ylim1,xaxt='n',col=col1,mgp=c(2.2,1,0))
  matpoints(c,wei,type="p",pch=2:6,main=main1,xlab=expression(c),ylab="Weight",lty=ltype1,lwd=lthick1,ylim=ylim1,xaxt='n',col=col1,mgp=c(2.2,1,0))
  axis(1,seq(0,1,0.2))
  legend('bottomright',c("RNA-op","RNA-2","RNA-5","RNA-10","RNA-n/2"), lty=ltype1,pch=2:6 ,col=col1,bg='white',cex=0.8)
}

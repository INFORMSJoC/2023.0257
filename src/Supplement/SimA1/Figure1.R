library(readxl)
SimA1_PL <- read_excel("C:/Users/dell/Desktop/2023.0257/scr/output/SimA1_PL.xlsx",col_names = FALSE)
View(SimA1_PL)

nme=7                                                                   
cx=c(0.01,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,0.99);                      
ns=c(100,200,500,1000)
nn=length(ns)
awei=array(0,dim=c(nn,length(cx),nme))

#Surrogate PL

#awei[1,,]=as.matrix(SimA1_PL[1:11,c(1:7)],11,7) 
#awei[2,,]=as.matrix(SimA1_PL[25:35,c(1:7)],11,7) 
#awei[3,,]=as.matrix(SimA1_PL[49:59,c(1:7)],11,7) 
#awei[4,,]=as.matrix(SimA1_PL[73:83,c(1:7)],11,7)  

#Targeted PL

awei[1,,]=as.matrix(SimA1_PL[12:22,c(1:7)],11,7) 
awei[2,,]=as.matrix(SimA1_PL[36:46,c(1:7)],11,7) 
awei[3,,]=as.matrix(SimA1_PL[60:70,c(1:7)],11,7) 
awei[4,,]=as.matrix(SimA1_PL[84:94,c(1:7)],11,7)  

lthick1=matrix(2,nme,1)
col1=c("green","red","purple","blue","black","brown","darkgreen")
ltype1=c(1,1,1,6,1,1,2)
ylim1=c(0.98,1.15)
#layout(matrix(c(1,2,3,4,5), ncol=2, byrow=TRUE), heights=c(4,4,0.5))
par(mfrow=c(2,2))
#pdf(file=figname,width = 6, height = 4 )
for (ni in 1:nn) {
  wei=awei[ni,,]
  main1 = bquote(paste('n','=',.(ns[ni])))
  par(mar=c(5, 5, 2, 1))
  matplot(cx,wei, type="l",main=main1,xlab=expression(c),ylab="Relative Targeted Empirical Ranking Risk",lty=ltype1,lwd=lthick1,ylim=ylim1,xaxt='n',col=col1,mgp=c(2.2,1,0))
  matpoints(cx,wei,type="p",pch=c(0,1,2,3,4,5,6),main=main1,xlab=expression(c),ylab="Relative Targeted Empirical Ranking Risk",lty=ltype1,lwd=lthick1,ylim=ylim1,xaxt='n',col=col1,mgp=c(2.2,1,0))
  axis(1,seq(0,1,0.2))
  #if(ni==1) 
  legend('topleft',c("Incorrect","Correct","RMA-op","RMA-2","RMA-5","RMA-10","RMA-n/2"), lty=ltype1,pch=c(0,1,2,3,4,5,6),col=col1,bg='white',cex=0.85)
}

#figname <- pase("figure_sample_",sam_size,".eps",sep="")
#savePlot(file=figname,type="eps",dev.cur())
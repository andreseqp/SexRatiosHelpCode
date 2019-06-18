# Figure A1 --------------------------------------------------------------------
# Replicates of the co-evolutionary dynamics of helping, 
# fratricide and sex ratios  


# Directory with the IBS -------------------------------------------------------
Simdir<-"E:\\Dropbox\\Haplodiploidy\\Gil\\Andres\\IBD"


png("E:\\Dropbox\\Haplodiploidy\\Gil\\Andres\\Evolution\\Figure_A1_dynamics.png",
    width=210*3,height=297*3)


parsize<-1.3
labsize<-2.5
# Colour ribbons
colorPol<-c("#a1d99b","#fed976","#fed976","paleturquoise3","#bcbddc")



list1<-list.files()
list1<-list1[grep('09',x=list1)]
listEv<-list1[grep('evol',x = list1)]
list.dist<-list1[grep('dist',x = list1)]


xaxs<-c("s",rep("n",4))
yaxs<-c("s","n")
ylabs<-cbind(c("","","phenotypic value","",""),rep("",5))
countx<-1
county<-6
plot.new()
for(seed in 1:10){
  if(county==1){county<-6;countx<-countx+1}
  county<-county-1
  Xs <- read.table(listEv[seed])
  names(Xs) <- c("Gen","z1","z2","z3","h","om","sdz1","sdz2","sdz3","sdh","sdom",
                 "Q1z1","Q1z2","Q1z3","Q1h","Q1om",
                 "Q3z1","Q3z2","Q3z3","Q3h","Q3om")
  par(las=1,plt=posPlot(numploty=5,idploty=county,numplotx=2,idplotx=countx),
      xpd=F,yaxt=yaxs[countx],xaxt=xaxs[county],new=T)
  plot(x=c(0,20000),y=c(0.5,0.5),type="l",lwd=2,col="grey",ylim=c(0,1),
       xlab="",ylab=ylabs[county,countx],xlim=c(0,15000),cex.lab=labsize,cex.axis=1.2)
  # Ribbons
  polygon(x=c(Xs$Gen,Xs$Gen[length(Xs$Gen):1]),
          y=c(Xs$Q3z1,Xs$Q1z1[length(Xs$sdz1):1]),
          col=colorPol[1],border=NA)
  polygon(x=c(Xs$Gen,Xs$Gen[length(Xs$Gen):1]),
          y=c(Xs$Q3z2,Xs$Q1z2[length(Xs$sdz2):1]),
          col=colorPol[2],border=NA)
  polygon(x=c(Xs$Gen,Xs$Gen[length(Xs$Gen):1]),
          y=c(Xs$Q3z3,Xs$Q1z3[length(Xs$sdz3):1]),
          col=colorPol[3],border=NA)
  polygon(x=c(Xs$Gen,Xs$Gen[length(Xs$Gen):1]),
          y=c(Xs$Q3h,Xs$Q1h[length(Xs$sdh):1]),
          col=colorPol[4],border=NA)
  polygon(x=c(Xs$Gen,Xs$Gen[length(Xs$Gen):1]),
          y=c(Xs$Q3om,Xs$om[length(Xs$om):1]-Xs$sdom[length(Xs$sdom):1]),
          col=colorPol[5],border=NA)
  lines(Xs$Gen,Xs$z1,lwd=3,col=colors1[1])
  lines(Xs$Gen,Xs$z2,lwd=3,col=colors1[2])
  lines(Xs$Gen,Xs$z3,lwd=3,col=colors1[3],lty=2)
  lines(Xs$Gen,Xs$h,lwd=3,col=colors1[4])
  lines(Xs$Gen,Xs$om,lwd=3,col=colors1[5])
  box()
  par(xpd=T)
}
text("Generations",x=-500,y=-0.565,cex=labsize)
legend("bottomright",legend=c("spring sex ratio","summer sex ratio solit.",
                           "summer sex ratio social.","helping","fratricide")
       ,lty=c(1,1,2,1,1),lwd=2,col=colors1,bty="o",cex=0.7)

dev.off()

# Parameters

text(x=10000,y=0.15,bquote(italic(B)==.(b*Sf),where=as.list(parameters)),cex=parsize)
text(x=14000,y=0.15,bquote(italic(rho)==.(phi),where=as.list(parameters)),cex=parsize)
text(y=0.95,x=500,"A",cex=labsize*0.7)

# Figure 4 ---------------------------------------------------------------------
# Minimal level of recycling efficiency required for the evolution of fratricide
# as a function of the sex ratio in social nests.

source("SelecDiff_Evol.R")

png(filename = "fig5_col.png",width = 765,height = 630)

# Plotting parameters 
lty1<-c(1,2)
colors1<-c("#e41a1c","#377eb8", "#4daf4a", "#984ea3")

# Parameters and plotting ranges
param<-c(F1=5,Sf=1)
sex.rat<-c(z1=0.5,z2=0.5)
SmRange<-seq(0,1,length=4)
HRange<-seq(0,1,length=4)
z3Range<-seq(0,1,length=50)
par(plt=posPlot(numplotx = 1,idplotx = 1),las=1,xaxt="s",yaxt="s",xpd=FALSE,
    fg="black",bg="white"
    ,xaxt="s",yaxt="s")
plot(y=c(1,1),x=c(0,1),xlim=c(0,1),ylim=c(0,1),col="white" ,lwd=1,ylab="",
     xlab="Sex-ratio in social nests",cex.axis=1.2,type="l",cex.lab=1.2,xaxt='n')
axis(1,at = c(0,0.3,0.5,0.7,1),cex.lab=1.4)
mtext(text=expression(italic(rho)[min]),2,line=3,cex=1.5)
for(i in 1:4){
  lines(Phim(SmRange[i],1,h=0.01,param,sex.rat,z3Range,0)~z3Range,
        col=colors1[i],lty=lty1[1],lwd=2)
}
lines(x=c(0.5,0.5),y=c(-0.1,1.1),col="grey")
legend("topleft",legend=c(expression(italic(S)[m]==0),
                          expression(italic(S)[m]==0.33),
                          expression(italic(S)[m]==0.66),
                          expression(italic(S)[m]==1)),
       col=colors1,lty=lty1[1],lwd=2,cex=1.3)

dev.off()


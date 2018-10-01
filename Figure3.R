# Figure 3 ---------------------------------------------------------------------
# Minimal worker benefits required for the emergence of helping behavior as 
# a function of fratricidal strategy (A) and the sex ratio (proportion sons) of 
# social nests (B). 

source("D:/quinonesa/Dropbox/Haplodiploidy/Gil/Andres/SexRatiosHelpCode/SelecDiff_Evol.R")

# Plotting parameters
lwd<-c(3,1.5)
colors1<-c("#e41a1c","#377eb8", "#4daf4a", "#984ea3")
colBck<-rgb(0,0,0,10,maxColorValue=255)
lty1<-c(1,2)

pngdir<-"D:\\Users\\quinonesa\\Dropbox\\Haplodiploidy\\Gil\\Andres\\Evolution\\"
png(paste(pngdir,"Figure_3_col.png"),width=1000,height=800)

# Panel A ----------------------------------------------------------------------

# Model parameters and initial conditions
param<-c(F1=5,Sf=1)
sex.rat<-c(z1=0.5,z2=0.5)
parameters <- c(F1=5,Sm=1,Sf=1,b=0,K=0.002,phi=0)
state <- c(z1=0.5,z2=0.5,z3=0.5,h=0,omega=0)
# Allow sex-ratio evolution
times = seq(0,100,by=1)
out = ode(y=state,times=times,func=Evolution,parms=parameters, method="lsode")
sex.rat.bias<-c(out[length(out[,1]),2],out[length(out[,1]),3])
OmRange<-seq(0,1,length=50)
phiRange<-seq(0,1,length=3)
par(plt=posPlot(numplotx = 2,idplotx = 1),las=1,xaxt="s",yaxt="s",xpd=FALSE,
    fg="black",bg="white"
    ,xaxt="s",yaxt="s")
plot(y=c(1,1),x=c(0,1),xlim=c(0,1),ylim=c(0.5,1.5),col="grey" ,lwd=1,ylab="",
     xlab="Fratricide",cex.axis=1.2,type="l",cex.lab=1.5)
par(las=1)
mtext(text=expression(italic(B)[min]),2,line=3,cex=1.5)
for(i in 1:4){
  lines(BminFuc(0,param,sex.rat,0.5,OmRange,phiRange[i])~OmRange,lwd=2,
        col=colors1[i],lty=lty1[1])
}
for(i in 1:4){
  lines(BminFuc(1,param,sex.rat.bias,out[length(out[,1]),3],OmRange,
                phiRange[i])~OmRange,lwd=2,col=colors1[i],lty=lty1[2])
}
leg<-legend("bottomright",legend=c(expression(paste(italic(rho)==0,"    ")),
                                   expression(italic(rho)==0.5),
                                   expression(italic(rho)==1),
                                   expression(paste(italic(rho)==0,"    ")),
                                   expression(italic(rho)==0.5),
                                   expression(italic(rho)==1))
            ,ncol=2,col=c(colors1[1:3],colors1[1:3]),lty=c(rep(lty1[1],3),
                                                           rep(lty1[2],3)),
            title=expression(paste(italic(S)[m]==0,"            ",
                                   italic(S)[m]==1)),bty="o",lwd=2,cex=1.2)
text(x=0.08,y=1.45,"A",cex=2)

# Panel B ----------------------------------------------------------------------

SmRange<-seq(0,1,length=4)
z3Range<-seq(0,1,length=50)
OmRange<-seq(0,1,length=4)
par(plt=posPlot(numplotx = 2,idplotx = 2),las=1,xaxt="s",yaxt="s",xpd=FALSE,
    fg="black",bg="white",xaxt="s",yaxt="n",new=T)
plot(y=c(1,1),x=c(0,1),xlim=c(0,1),ylim=c(0.5,1.5),col="grey" ,lwd=1,ylab="",
     xlab="Sex ratio in social nests",cex.axis=1.2,type="l",cex.lab=1.5)
par(las=1)
for(i in 1:4){
  lines(BminFuc(SmRange[i],param,sex.rat,z3Range,0,0.6)~z3Range,lwd=2,
        col=colors1[i],lty=lty1[1])
}
legend("bottomright",legend=c(expression(italic(S)[m]==0),
                              expression(italic(S)[m]==1/3),
                              expression(italic(S)[m]==2/3),
                              expression(italic(S)[m]==1)),
       col=colors1,lty=lty1[1],lwd=2,cex=1.2)
text(x=0.08,y=1.45,"B",cex=2)

dev.off()

# Figure 5 ---------------------------------------------------------------------
# Strength and direction of selection in the summer sex ratios of social (A) 
# and solitary nests (B)

source("D:/quinonesa/Dropbox/Haplodiploidy/Gil/Andres/SexRatiosHelpCode/SelecDiff_Evol.R")

png(filename = paste(pngdir,Sys.Date(),"fig6_col.png",sep=""),width = 765,
    height = 630)

# Plotting parameters
colors1<-c("#e41a1c","#377eb8", "#4daf4a", "#984ea3")
lty1<-c(1,2)

# Panel A ----------------------------------------------------------------------

Bplot<-1
h<-0.1
phirange<-c(0,0.3,0.6,0.9)
omegRange<-seq(0,1,length=50)
param<-c(F1=5,Sf=1,Sm=1)
param2<-c(F1=5,Sf=1,Sm=1)
par(plt=posPlot(numplotx=2,idplotx=1)+c(-0.047,-0.052,0,0),las=1,xaxt="s",
    yaxt="s",xpd=FALSE,fg="black",bg="white",xaxt="s")
plot(x=c(0,2),y=c(0,0),ylim=c(-0.2,0.1),xlim=c(0,1),col="grey" ,lwd=1,ylab="",
     xlab="Fratricide",cex.axis=1,type="l",cex.lab=1.2)
for(i in 1:4){
  lines(diffz3(Bplot,param2,z1=0.5,z2=0.5,z3=0.5,omegRange,h,phirange[i])~
          omegRange,col=colors1[i],lty=lty1[1],lwd=2)
}
par(xpd=T,las=0)
text(-0.28,-0.05,expression(paste("Selection on ",italic(z[3]))),srt=90,cex=1.2)

parameters <- c(F1=5,Sm=1,Sf=1,b=0,K=0.002,phi=0)
state <- c(z1=0.5,z2=0.5,z3=0.5,h=0,omega=0)
# Create output
times = seq(0,10000,by=10)
out = ode(y=state,times=times,func=Evolution,parms=parameters, method="lsode")

for(j in 1:4){
  lines(diffz3(Bplot,param2,z1=out[1001,2],z2=out[1001,3],z3=out[1001,3],
               omegRange,h,phirange[j])~omegRange,
        col=colors1[j],lty=lty1[2],lwd=2)
}

legend("topright",legend=c(expression(paste(rho==0,"       ")),
                           expression(rho==0.3),
                           expression(rho==0.6),
                           expression(rho==0.9),
                           expression(paste(rho==0,"       ")),
                           expression(rho==0.3),
                           expression(rho==0.6),
                           expression(rho==0.9)),
       ncol=2,lty=c(rep(lty1[1],4),rep(lty1[2],4)),lwd=2,col=colors1,bty="o",
       title=paste("Even sex-ratios","  ","Equilibrium sex-ratios"))
par(bg="white")
text(x=0.05,y=0.085,"A",cex=2)

# Panel B ----------------------------------------------------------------------

Bplot<-1
hrange<-seq(0,0.3,length=3)
z3Range<-seq(0,1,length=100)
phiplot<-0.8
omegaPlot<-0.8
param<-c(F1=5,Sf=1,Sm=1)
param2<-c(F1=5,Sf=1,Sm=1,phi=phiplot)
par(plt=posPlot(numplotx=2,idplotx=2)+c(-0.052,-0.047,0,0),las=1,xaxt="s",
    yaxt="s",xpd=FALSE,fg="black",bg="white",xaxt="s",new=T,cex.axis=1)
plot(x=c(-0.2,1.1),y=c(0,0),ylim=c(-1,1),xlim=c(0,1),col="grey" ,lwd=1,ylab="",
     xlab="Sex-ratio in social nests",
     cex.axis=1,type="l",cex.lab=1.2,yaxt="n")
axis(side=4)
lines(x=c(0.5,0.5),y=c(-1.1,1.1),col="grey")
for(i in 1:3){
  lines(diffz2(Bplot,param2,z1=0.5,z2=0.5,z3=z3Range,omegaPlot,hrange[i])~
          z3Range,col=colors1[i],lty=lty1[1],lwd=2)
}

par(xpd=T,las=0)
text(1.2,0.0,expression(paste("Selection on ",italic(z[2]))),srt=90,cex=1.2)
parameters <- c(F1=5,Sm=1,Sf=1,b=0,K=0.002,phi=0)
state <- c(z1=0.5,z2=0.5,z3=0.5,h=0,omega=0)
# Create output
times = seq(0,10000,by=10)
out = ode(y=state,times=times,func=Evolution,parms=parameters, method="lsode")

for(j in 1:3){
  lines(diffz2(Bplot,param2,z1=out[1001,2],z2=out[1001,3],z3=z3Range,omegaPlot,
               hrange[j])~z3Range,col=colors1[j],lty=lty1[2],lwd=2)
}

# par(bg=colBck)
legend("topright",legend=c(expression(paste(italic(h)==0,"       ")),
                           expression(italic(h)==0.1),
                           expression(italic(h)==0.2),
                           expression(paste(italic(h)==0,"       ")),
                           expression(italic(h)==0.1),
                           expression(italic(h)==0.2)),
       ncol=2,lty=c(rep(lty1[1],3),rep(lty1[2],3)),lwd=2,col=c(colors1[1:3],
                                                               colors1[1:3]),
       bty="o",title=paste("Even sex-ratios","  ","Equilibrium sex-ratios"))
par(bg="white")
text(x=0.05,y=0.9,"B",cex=2)

dev.off()

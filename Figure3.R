# Figure 3 ---------------------------------------------------------------------
# Minimal worker benefits required for the emergence of helping behavior as 
# a function of fratricidal strategy (A) and the sex ratio (proportion sons) of 
# social nests (B). 

source("SelecDiff_Evol.R")
library("data.table")

# Plotting parameters
lwd<-c(3,1.5)
colors1<-c("#e41a1c","#377eb8", "#4daf4a", "#984ea3")
colBck<-rgb(0,0,0,10,maxColorValue=255)
lty1<-c(1,2)


png("Figure_3_col.png",width=1000,height=800)

# Panel A ----------------------------------------------------------------------

# Model parameters and initial conditions
lty1<-c(1,2)
param<-c(F1=5,Sf=1)
sex.rat<-c(z1=0.5,z3=0.5)
parameters <- c(F1=5,Sm=1,Sf=1,b=0,K=0.002,phi=0)
state <- c(z1=0.5,z3=0.5,z5=0.5,h=0,omega=0)
# Create output
times = seq(0,100,by=1)
out = ode(y=state,times=times,func=Evolution,parms=parameters, method="lsode")
sex.rat.bias<-c(out[length(out[,1]),2],out[length(out[,1]),3])
OmRange<-seq(0,1,length=50)
phiRange<-c(0,1)
z5Range<-c(0,0.2,0.5,0.75)

par(plt=posPlot(numplotx = 2,idplotx = 1)+c(-0.047,-0.052,0,0),
    las=1,xaxt="s",yaxt="s",xpd=FALSE,fg="black",bg="white"
    ,xaxt="s",yaxt="s")
plot(y=c(1,1),x=c(0,1),xlim=c(0,1),ylim=c(0.3,2),col="grey" ,lwd=1,ylab="",
     xlab="Fratricide",cex.axis=1.2,type="l",cex.lab=1.5)
par(las=1)
mtext(text=expression(italic(B)[min]),2,line=3,cex=1.5)
for(i in 1:4){
  lines(BminFuc(1,param,sex.rat.bias,z5Range[i],OmRange,phiRange[1])~OmRange,lwd=2,
        col=colors1[i],lty=lty1[1])
}
for(i in 1:4){
  lines(BminFuc(1,param,sex.rat.bias,z5Range[i],OmRange,phiRange[2])~OmRange,
        lwd=2,col=colors1[i],lty=lty1[2])
}

points(x=0,y=0.72,pch=20,cex=1.5)
points(x=0.5,y=0.72,pch=20,cex=1.5)
arrows(x0 = 0,x1 = 0.5,y0 = 0.72,y1 = 0.72,lwd = 1.2,length = 0.1)
# x=c(-0.05,0.15),y=c(1.4,1.5).x=c(0.2,0.3),y=c(1.4,1.5)
leg<-legend("bottomright",legend=c(expression(paste(italic(z[5])==0,"    ")),
                                   expression(italic(z[5]==0.2)),
                                   expression(italic(z[5])==0.5),
                                   expression(italic(z[5])==0.75),
                                   expression(paste(italic(z[5])==0,"    ")),
                                   expression(italic(z[5])==0.2),
                                   expression(italic(z[5])==0.5),
                                   expression(italic(z[5])==0.75))
            ,ncol=2,col=c(colors1[1:4],colors1[1:4]),lty=c(rep(lty1[1],4),
                                                           rep(lty1[2],4)),
            title=expression(paste(italic(rho)==0,"                     ",
                                   italic(rho)==1)),bty="o",lwd=2,cex=1.2)
text(x=par('usr')[1]+(par('usr')[2]-par('usr')[1])*0.1,
     y=par('usr')[4]-(par('usr')[4]-par('usr')[3])*0.1,"A",cex=2)
# Panel B ----------------------------------------------------------------------

# second panel. Long term effect of fratricide in helping behaviour

# Create output
parameters <- c(F1=5,Sm=1,Sf=1,b=0.85,K=0.005,phi=0.7)
state <- c(z1=0.5,z3=0.5,z5=0.5,h=0,omega=0.0)
times = seq(0,20000,by=10)
eventdat = data.frame(var=c("h","omega"),time=c(1000,1000),
                      value=c(0.0002,0.0002),method=c("add"))
eventdat2 = data.frame(var=c("h"),time=c(1000),
                       value=c(0.0002),method=c("add"))

bRange<-seq(0.7,1.2,length=100)
rhoRange<-c(0.6,0.7,0.8)

# Numerical simulations run on fresh start

LastGenFrat<-data.table(matrix(data=0,nrow = 300,ncol = 7))
LastGen<-data.table(matrix(data=0,nrow = 300,ncol = 7))

names(LastGen)<-c("b",'rho','z1','z3','z5','h','ome')
names(LastGenFrat)<-c("b",'rho','z1','z3','z5','h','ome')

LastGenFrat[,'b']<-rep(bRange,3)
LastGenFrat[,'rho']<-rep(rhoRange,each=100)
LastGen[,'b']<-rep(bRange,3)
LastGen[,'rho']<-rep(rhoRange,each=100)

for(rhoit in rhoRange){
  parameters["phi"]<-rhoit
  for (bit in bRange) {
    parameters["b"]<-bit
    out <-ode(y=state,times=times,func=Evolution,
              parms=parameters,events=list(data=eventdat), method="lsode")

    out2 <-ode(y=state,times=times,func=Evolution,
               parms=parameters,events=list(data=eventdat2), method="lsode")
    LastGenFrat[b==bit&rho==rhoit,':='(z1=out[dim(out)[1],2],
                                       z3=out[dim(out)[1],3],
                                       z5=out[dim(out)[1],4],
                                       h=out[dim(out)[1],5],
                                       ome=out[dim(out)[1],6])]
    LastGen[b==bit&rho==rhoit,':='(z1=out2[dim(out2)[1],2],
                                   z3=out2[dim(out2)[1],3],
                                   z5=out2[dim(out2)[1],4],
                                   h=out2[dim(out2)[1],5],
                                   ome=out2[dim(out2)[1],6])]
  }
}

par(las=1,plt=posPlot(numploty=1,numplotx = 2,idplotx = 2,idploty=1)
    +c(-0.047,-0.052,0,0),xpd=F,new=T)
plot(h~b,data = LastGen[rho==0.6],type="l",lwd=3,col=colors1[1],ylim=c(0,1),
     xlim=c(0.7,1.2),xlab="benefits of helping",ylab="",xaxt='s',yaxt='n',
     cex.lab=1.5,cex.axis=1.2,lty=lty1[1])
axis(4)
mtext(text = "Evolved helping level",side = 4,las=3,line = 3,cex=2)
lines(h~b,data = LastGenFrat[rho==0.7],lwd=3,col=colors1[2],lty=lty1[1])
lines(x = c(1,1),y=c(0,1),col='grey')
legend("bottomright",legend=c('without fratricide','with fratricide')
       ,lty=c(1,1,2,1,1),lwd=2,col=colors1[1:2],bty="o",cex=1.2)
text(x=par('usr')[1]+(par('usr')[2]-par('usr')[1])*0.1,
     y=par('usr')[4]-(par('usr')[4]-par('usr')[3])*0.1,"B",cex=2)


dev.off()

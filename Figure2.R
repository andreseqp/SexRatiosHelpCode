# Figure 2 --------------------------------------------------------------------
# Co-evolutionary dynamics of helping, fratricide and sex ratios. 

library(here)
source("SelecDiff_Evol.R")

# Directory with the IBS -------------------------------------------------------


png("Figure_2_dynamics.png",width=1000,height=800)


  parsize<-1.3
  labsize<-2.5
  
  # Panel B ----------------------------------------------------------------------
  
  pertTime<-5000
  
  parameters <- c(F1=5,Sm=0.8,Sf=0.8,b=0.9,K=0.005,phi=0.7)
  #Intitial conditions
  state <- c(z1=0.5,z3=0.5,z5=0.5,h=0,omega=0.0)
  # Create output of numerical simulations 
  times = seq(0,20000,by=10)
  # Let helping and fratricide start evolving 
  eventdat = data.frame(var=c("h","omega"),time=c(1000,1000),
                        value=c(0.0002,0.0002),method=c("add"))
  # Perturb equilibria in the direction of IBS
  eventdat2 = data.frame(var=c("h","omega","z3","z5"),time=rep(pertTime,4),
                         value=c(0.0002,0.502,0.0,-0.2),method=c("add"))
  out = ode(y=state,times=times,func=Evolution,parms=parameters,
            events=list(data=eventdat2), method="lsode")
  
  
  colors1<-c("#31a354","#d95f0e","#d95f0e","#3182bd","Purple")
  colBck<-rgb(0,0,0,10,maxColorValue=255)
  
  
  
  list1<-list.files(here('IBD'))
  list1<-list1[grep('09',x=list1)]
  listEv<-list1[grep('evol',x = list1)]
  list.dist<-list1[grep('dist',x = list1)]
  seed<-2
  
  Xs <- read.table(here('IBD',listEv[seed]))
  names(Xs) <- c("Gen","z1","z3","z5","h","om","sdz1","sdz3","sdz5","sdh","sdom",
                 "Q1z1","Q1z3","Q1z5","Q1h","Q1om",
                 "Q3z1","Q3z3","Q3z5","Q3h","Q3om")
  
  
  par(las=1,plt=posPlot(numploty=2,numplotx=2,idploty=1,idplotx=1),xpd=F)
  plot(out[,1],out[,2],type="l",lwd=3,col=colors1[1],ylim=c(0,1),xlim=c(0,15000),
       xlab="",ylab="",cex.lab=1.5,cex.axis=1.2,xaxt="s",yaxt="s")
  abline(h=0.5,lwd=2,col="grey")
  lines(out[,1],out[,3],lwd=3,col=colors1[2])
  lines(out[,1],out[,4],lwd=3,col=colors1[3],lty=2)
  lines(out[,1],out[,5],lwd=3,col=colors1[4])
  lines(out[,1],out[,6],lwd=3,col=colors1[5])
  par(bg="white",xpd=TRUE)
  
  text(x=-4500,y=1.2,labels="phenotypic value",cex=labsize,srt=90)
  
  par(bg="white",xpd=FALSE)
  text(x=10000,y=0.15,bquote(italic(B)==.(b*Sf),where=as.list(parameters)),
       cex=parsize)
  text(x=14000,y=0.15,bquote(italic(rho)==.(phi),where=as.list(parameters)),
       cex=parsize)
  text(y=0.95,x=500,"B",cex=labsize*0.7)
  
  # Panel A ----------------------------------------------------------------------
  
  # Colour ribbons
  colorPol<-c("#a1d99b","#fed976","#fed976","paleturquoise3","#bcbddc")
  par(las=1,plt=posPlot(numploty=2,idploty=2,numplotx=2,idplotx=1),xpd=F,
      new=T,yaxt="s",xaxt='n')
  plot(x=c(0,20000),y=c(0.5,0.5),type="l",lwd=2,col="grey",ylim=c(0,1),
       xlab="",ylab="",xlim=c(0,15000),cex.lab=1.5,cex.axis=1.2)
  # Ribbons
  polygon(x=c(Xs$Gen,Xs$Gen[length(Xs$Gen):1]),
          y=c(Xs$Q3z1,Xs$Q1z1[length(Xs$sdz1):1]),
          col=colorPol[1],border=NA)
  polygon(x=c(Xs$Gen,Xs$Gen[length(Xs$Gen):1]),
          y=c(Xs$Q3z3,Xs$Q1z3[length(Xs$sdz3):1]),
          col=colorPol[2],border=NA)
  polygon(x=c(Xs$Gen,Xs$Gen[length(Xs$Gen):1]),
          y=c(Xs$Q3z5,Xs$Q1z5[length(Xs$sdz5):1]),
          col=colorPol[3],border=NA)
  polygon(x=c(Xs$Gen,Xs$Gen[length(Xs$Gen):1]),
          y=c(Xs$Q3h,Xs$Q1h[length(Xs$sdh):1]),
          col=colorPol[4],border=NA)
  polygon(x=c(Xs$Gen,Xs$Gen[length(Xs$Gen):1]),
          y=c(Xs$Q3om,Xs$om[length(Xs$om):1]-Xs$sdom[length(Xs$sdom):1]),
          col=colorPol[5],border=NA)
  
  # Numerical simulations
  lines(Xs$Gen,Xs$z1,lwd=3,col=colors1[1])
  lines(Xs$Gen,Xs$z3,lwd=3,col=colors1[2])
  lines(Xs$Gen,Xs$z5,lwd=3,col=colors1[3],lty=2)
  lines(Xs$Gen,Xs$h,lwd=3,col=colors1[4])
  lines(Xs$Gen,Xs$om,lwd=3,col=colors1[5])
  box()
  par(xpd=T)
  
  text(x=10000,y=0.15,bquote(italic(B)==.(b*Sf),where=as.list(parameters)),cex=parsize)
  text(x=14000,y=0.15,bquote(italic(rho)==.(phi),where=as.list(parameters)),cex=parsize)
  text(y=0.95,x=500,"A",cex=labsize*0.7)
  
  # Panel C ----------------------------------------------------------------------
  
  parameters <- c(F1=5,Sm=0.8,Sf=0.8,b=1.2,K=0.005,phi=0.7)
  # Initial conditions
  state <- c(z1=0.5,z3=0.5,z5=0.5,h=0,omega=0.0)
  # Create output
  times = seq(0,20000,by=10)
  eventdat = data.frame(var=c("h","omega"),time=c(1000,1000),
                        value=c(0.0002,0.0002),method=c("add"))
  out = ode(y=state,times=times,func=Evolution,parms=parameters,
            events=list(data=eventdat), method="lsode")
  
  
  colors1<-c("#31a354","#d95f0e","#d95f0e","#3182bd","Purple")
  colBck<-rgb(0,0,0,10,maxColorValue=255)
  
  par(las=1,plt=posPlot(numploty=2,idploty=2,numplotx=2,idplotx=2),xpd=F,new=T)
  plot(out[,1],out[,2],type="l",lwd=3,col=colors1[1],ylim=c(0,1),xlim=c(0,8000),
       xlab="",ylab="",cex.lab=1.5,cex.axis=1.2,xaxt="n",yaxt="n")
  abline(h=0.5,lwd=2,col="grey")
  lines(out[,1],out[,3],lwd=3,col=colors1[2])
  lines(out[,1],out[,4],lwd=3,col=colors1[3],lty=2)
  lines(out[,1],out[,5],lwd=3,col=colors1[4])
  lines(out[,1],out[,6],lwd=3,col=colors1[5])
  par(bg=colBck)
  par(bg="white")
  
  legend("topright",legend=c("spring sex ratio","summer sex ratio solit.",
                             "summer sex ratio social.","helping","fratricide")
         ,lty=c(1,1,2,1,1),lwd=2,col=colors1,bty="o",cex=1.15)
  
  text(x=4000,y=0.15,bquote(italic(B)==.(b*Sf),where=as.list(parameters)),
       cex=parsize)
  text(x=6500,y=0.15,bquote(italic(rho)==.(phi),where=as.list(parameters)),
       cex=parsize)
  text(y=0.95,x=250,"C",cex=labsize*0.7)
  
  
  # Panel D ----------------------------------------------------------------------
  
  list1<-list.files(here('IBD'))
  list1<-list1[grep('15',x=list1)]
  listEv<-list1[grep('evol',x = list1)]
  list.dist<-list1[grep('dist',x = list1)]
  seed<-3
  Xs <- read.table(here('IBD',listEv[seed]))
  names(Xs) <- c("Gen","z1","z3","z5","h","om","sdz1","sdz3","sdz5","sdh","sdom",
                 "Q1z1","Q1z3","Q1z5","Q1h","Q1om",
                 "Q3z1","Q3z3","Q3z5","Q3h","Q3om")
  
  dist<-read.table(here('IBD',list.dist[seed]))
  names(dist) <- c("z1","z3","z5","h","om")
  
  colorPol<-c("#a1d99b","#fed976","#fed976","paleturquoise3","#bcbddc")
  par(las=1,plt=posPlot(numploty=2,idploty=1,numplotx=2,idplotx=2),xpd=F,new=T,yaxt="n")
  plot(x=c(0,20000),y=c(0.5,0.5),type="l",lwd=2,col="grey",ylim=c(0,1),
       xlab="",ylab="",xlim=c(0,8000),cex.lab=1.5,cex.axis=1.2,xaxt='s')
  #  SD
  # polygon(x=c(Xs$Gen,Xs$Gen[length(Xs$Gen):1]),y=c(Xs$z3+Xs$sdz3,Xs$z3[length(Xs$z3):1]-Xs$sdz3[length(Xs$sdz3):1]),
  #         col=colorPol[2],border=NA)
  # polygon(x=c(Xs$Gen,Xs$Gen[length(Xs$Gen):1]),y=c(Xs$z5+Xs$sdz5,Xs$z5[length(Xs$z1):1]-Xs$sdz5[length(Xs$sdz5):1]),
  #         col=colorPol[3],border=NA)
  # polygon(x=c(Xs$Gen,Xs$Gen[length(Xs$Gen):1]),y=c(Xs$h+Xs$sdh,Xs$h[length(Xs$h):1]-Xs$sdh[length(Xs$sdh):1]),
  #         col=colorPol[4],border=NA)
  # polygon(x=c(Xs$Gen,Xs$Gen[length(Xs$Gen):1]),y=c(Xs$om+Xs$sdom,Xs$om[length(Xs$om):1]-Xs$sdom[length(Xs$sdom):1]),
  #         col=colorPol[5],border=NA)
  
  # IQR
  polygon(x=c(Xs$Gen,Xs$Gen[length(Xs$Gen):1]),
          y=c(Xs$z1+Xs$sdz1,Xs$z1[length(Xs$z1):1]-Xs$sdz1[length(Xs$sdz1):1]),
          col=colorPol[1],border=NA)
  polygon(x=c(Xs$Gen,Xs$Gen[length(Xs$Gen):1]),
          y=c(Xs$Q3z5,Xs$Q1z5[length(Xs$sdz5):1]),
          col=colorPol[3],border=NA)
  polygon(x=c(Xs$Gen,Xs$Gen[length(Xs$Gen):1]),
          y=c(Xs$Q3h,Xs$Q1h[length(Xs$sdh):1]),
          col=colorPol[4],border=NA)
  polygon(x=c(Xs$Gen,Xs$Gen[length(Xs$Gen):1]),
          y=c(Xs$Q3om,Xs$om[length(Xs$om):1]-Xs$sdom[length(Xs$sdom):1]),
          col=colorPol[5],border=NA)
  lines(Xs$Gen,Xs$z1,lwd=3,col=colors1[1])
  # lines(Xs$Gen,Xs$z3,lwd=3,col=colors1[2])
  lines(Xs$Gen,Xs$z5,lwd=3,col=colors1[3],lty=2)
  lines(Xs$Gen,Xs$h,lwd=3,col=colors1[4])
  lines(Xs$Gen,Xs$om,lwd=3,col=colors1[5])
  lines(x=c(0,20000),y=c(0.5,0.5),lwd=2,col="grey")
  box()
  par(xpd=T)
  text("Generations",x=-500,y=-0.265,cex=labsize)
  text(y=0.95,x=250,"D",cex=labsize*0.7)
  
  text(x=4000,y=0.15,expression(italic(B)==1.2),cex=parsize)
  text(x=6500,y=0.15,expression(italic(rho)==0.7),cex=parsize)
  
  # Inset panel D ----------------------------------------------------------------
  
  par(plt=posPlot(numplotx=6,numploty=6,idplotx=6,idploty=2)-0.01,new=T,xpd=T)
  hist(dist$z5,breaks=30,col=colorPol[3],main="",ylab="",xlab="",axes=F)
  polygon(x=c(-0.05,1.05,1.05,-0.05,-0.05),y=c(-10,-10,1000,1000,-10), 
          col = "white",border=T)
  par(plt=posPlot(numplotx=6,numploty=6,idplotx=6,idploty=2)-0.01,new=T)
  hist(dist$z5,breaks=30,col=colorPol[3],main="",ylab="",xlab="",axes=F)
  text(x=0.75,y=700,expression(z[5]),cex=2)

dev.off()

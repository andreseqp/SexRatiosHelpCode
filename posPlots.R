posPlot<-function(numplotx=1,numploty=1,idplotx=1,idploty=1,title=F,lowboundx=15)
{
  if(idplotx>numplotx){print("error x id does not match number of plots")}
  else{
  if(idploty>numploty){print("error y id does not match number of plots")}
  else{
  upboundx<-98
  #lowboundx<-15
  if(title){upboundy<-95}
  else {upboundy<-98}
  lowboundy<-15
  upxs<--1#rep(-1,2)numplotx)
  lowxs<-1#rep(1,2)numplotx)
  upys<--1#rep(-1,2)numploty)
  lowys<-1#rep(1,2)numploty)
  stepx<-(upboundx-lowboundx)/numplotx
  stepy<-(upboundy-lowboundy)/numploty
  upxs<-lowboundx+idplotx*stepx+upxs*0.5
  lowxs<-lowboundx+(idplotx-1)*stepx+lowxs*0.5
  upys<-lowboundy+idploty*stepy+upys*0.5
  lowys<-lowboundy+(idploty-1)*stepy+lowys*0.5
  return(c(lowxs,upxs,lowys,upys)/100)}}
}

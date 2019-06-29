# Functions necessary to produce figures ---------------------------------------

# Selection differentials ------------------------------------------------------

library(deSolve)
source('posPlots.R')

diffz1<-function(b,param,z1,z3,z5,omega,h)
{
  with(as.list(c(param)),
       {
         ht<-1-h
         z1t<-1-z1
         z3t<-1-z3
         z5t<-1-z5
         F3<-F1
         omegat<-1-omega
         H<-h*z1t*F1
         pi<-1-exp(-H)
         pit<-1-pi
         F5<-F3*(1+b*H)
         B<-Sf*b
         u3<-ht*z1t*F1
         u5<-Sf
         
         ZA<-pi*z5*(omegat+z5*omega*phi)+pit*z3
         
         ZAt<-pi*z5t*(1+z5*omega*phi)+pit*z3t
         
         V<-(1/2)*ht*z3t*F3
         
         G<-Sf*F3*pit+(1/2)*ht*z1t*F1*F3
         
         L<-Sf*F5*pi*z5*(omegat+z5*omega*phi)
         
         Om<-(Sm*z1*F1)/(z3*G+Sm*z1*F1+L)
         
         Im<-Sf*F5*pi*((1+z5*phi*omega)*((z3-z5)/z3t)+z5*omega)/(z3*G+Sm*z1*F1+L)
         
         
         diff<-F1*(V*(-1+(z1t/z1)+(1/2)*Om-(1/2)*Im)+(1/2)*Sm*(z3t/z3)*(1-Om+Im)
                   -(1/2)*Sf*h*F3*((pi*b+pit*(1+b*F1*h*z1t))*
                                     (z5t*(1+z5*omega*phi)+(z3t/z3)*z5*(omegat+z5*omega*phi)*(1-Om+Im))
                                   -pit*z3t*(2-Om+Im)))
         return(diff)
       })
}


diffz3<-function(b,param,z1,z3,z5,omega,h)
{
  with(as.list(c(param)),
       {
         
         ht<-1-h
         z1t<-1-z1
         z3t<-1-z3
         z5t<-1-z5
         F3<-F1
         omegat<-1-omega
         H<-h*z1t*F1
         pi<-1-exp(-H)
         pit<-1-pi
         F5<-F3*(1+b*H)
         B<-Sf*b
         u3<-ht*z1t*F1
         u5<-Sf
         
         ZA<-pi*z5*(omegat+z5*omega*phi)+pit*z3
         
         ZAt<-pi*z5t*(1+z5*omega*phi)+pit*z3t
         
         V<-(1/2)*ht*z3t*F3
         
         G<-Sf*F3*pit+(1/2)*ht*z1t*F1*F3
         
         L<-Sf*F5*pi*z5*(omegat+z5*omega*phi)
         
         Om<-(Sm*z1*F1)/(z3*G+Sm*z1*F1+L)
         
         Im<-Sf*F5*pi*((1+z5*phi*omega)*((z3-z5)/z3t)+z5*omega)/(z3*G+Sm*z1*F1+L)
         
         (vf1<-V*(2-Om+Im))
         (vm1<-V*(z1t/z1)+(1/2)*Sm*(z3t/z3)*(1-Om+Im))
         (vm2<-(1/2)*(z3t/z3)*(1-Om+Im))
         
         ## Selection differentials
         
         diff<-(1/2)*((u3*F3+u5*F5*pit)/(u3+u5))*(-1+(z3t/z3)*(1-Om+Im))
         
         return(diff)
       })
}


diffz5<-function(b,param,z1,z3,z5,omega,h,phi)
{
  with(as.list(c(param)),
       {
         
         ht<-1-h
         z1t<-1-z1
         z3t<-1-z3
         z5t<-1-z5
         F3<-F1
         omegat<-1-omega
         H<-h*z1t*F1
         pi<-1-exp(-H)
         pit<-1-pi
         F5<-F3*(1+b*H)
         B<-Sf*b
         u3<-ht*z1t*F1
         u5<-Sf
         
         ZA<-pi*z5*(omegat+z5*omega*phi)+pit*z3
         
         ZAt<-pi*z5t*(1+z5*omega*phi)+pit*z3t
         
         V<-(1/2)*ht*z3t*F3
         
         G<-Sf*F3*pit+(1/2)*ht*z1t*F1*F3
         
         L<-Sf*F5*pi*z5*(omegat+z5*omega*phi)
         
         Om<-(Sm*z1*F1)/(z3*G+Sm*z1*F1+L)
         
         Im<-Sf*F5*pi*((1+z5*phi*omega)*((z3-z5)/z3t)+z5*omega)/(z3*G+Sm*z1*F1+L)
         
         #   (vf1<-V*(2-Om+Im))
         #   (vm1<-V*(z1t/z1)+(1/2)*Sm*(z3t/z3)*(1-Om+Im))
         #   (vm2<-(1/2)*(z3t/z3)*(1-Om+Im))
         
         ## Selection differentials
         
         diff<-((1/2)*Sf*F5*pi/(u3+u5))*(-1+omega*phi*(1-2*z5)+(z3t/z3)*(1-Om+Im)*(omegat+2*omega*phi*z5))  
         
         return(diff)
       })
}



diffH<-function(b,param,z1,z3,z5,omega,h,phi)
{
  with(as.list(c(param)),
       {
         
         ht<-1-h
         z1t<-1-z1
         z3t<-1-z3
         z5t<-1-z5
         F3<-F1
         omegat<-1-omega
         H<-h*z1t*F1
         pi<-1-exp(-H)
         pit<-1-pi
         F5<-F3*(1+b*H)
         B<-Sf*b
         
         G<-Sf*F3*pit+(1/2)*ht*z1t*F1*F3
         
         L<-Sf*F5*pi*z5*(omegat+z5*omega*phi)
         
         Om<-(Sm*z1*F1)/(z3*G+Sm*z1*F1+L)
         
         Im<-Sf*F5*pi*((1+z5*phi*omega)*((z3-z5)/z3t)+z5*omega)/(z3*G+Sm*z1*F1+L)
         
         ## Selection differentials
         
         diff<-(1/2)*F5*(-z3t*(2-Om+Im)+(1/2)*B*(3*z5t*(1+z5*omega*phi)
                                                 +z5*(omegat+z5*omega*phi)*(z3t/z3)*(1-Om+Im)))
         
         return(diff)
       })
}

BminFuc<-function(Sm,param,sex.rat,z5,omega,phi){
  with(as.list(c(param,sex.rat)),
       {
         h<-0
         ht<-1-h
         z1t<-1-z1
         z3t<-1-z3
         z5t<-1-z5
         F3<-F1
         omegat<-1-omega
         H<-h*z1t*F1
         pi<-1-exp(-H)
         pit<-1-pi
         F5<-F3
         G<-Sf*F3*pit+(1/2)*ht*z1t*F1*F3
         L<-Sf*F5*pi*z5*(omegat+z5*omega*phi)
         Om<-(Sm*z1*F1)/(z3*G+Sm*z1*F1+L)
         Im<-Sf*F5*pi*((1+z5*phi*omega)*((z3-z5)/z3t)+z5*omega)/(z3*G+Sm*z1*F1+L)
         Bmin<-(2*z3t*(2-Om+Im))/
           (3*z5t*(1+z5*omega*phi)+z5*(omegat+z5*omega*phi)*(z3t/z3)*(1-Om+Im))
         return(Bmin)
       })
}

diffOm<-function(b,param,z1,z3,z5,h,phi)
{
  with(as.list(c(param)),
       {
         
         ht<-1-h
         z1t<-1-z1
         z3t<-1-z3
         z5t<-1-z5
         F3<-F1
         omegat<-1-omega
         H<-h*z1t*F1
         pi<-1-exp(-H)
         pit<-1-pi
         F5<-F3*(1+b*H)
         B<-Sf*b
         u3<-ht*z1t*F1
         u5<-Sf
         
         ZA<-pi*z5*(omegat+z5*omega*phi)+pit*z3
         
         ZAt<-pi*z5t*(1+z5*omega*phi)+pit*z3t
         
         V<-(1/2)*ht*z3t*F3
         
         G<-Sf*F3*pit+(1/2)*ht*z1t*F1*F3
         
         L<-Sf*F5*pi*z5*(omegat+z5*omega*phi)
         
         Om<-(Sm*z1*F1)/(z3*G+Sm*z1*F1+L)
         
         Im<-Sf*F5*pi*((1+z5*phi*omega)*((z3-z5)/z3t)+z5*omega)/(z3*G+Sm*z1*F1+L)
         
         #   (vf1<-V*(2-Om+Im))
         #   (vm1<-V*(z1t/z1)+(1/2)*Sm*(z3t/z3)*(1-Om+Im))
         #   (vm2<-(1/2)*(z3t/z3)*(1-Om+Im))
         
         ## Selection differentials
         
         #(omediff<-K*(1/4)*Sf*F5*h*B*z5*(phi*(3*z5t+z5*(z3t/z3)*(1-Om+Im))-z5*(z3t/z3)*(1-Om+Im)))/(u3+u5)
         diff<-(1/4)*Sf*F5*h*B*z5*(phi*(3*z5t+z5*(z3t/z3)*(1-Om+Im))-(z3t/z3)*(1-Om+Im))/(u3+u5)
         
         
         
         return(diff)
       })
}

Phim<-function(Sm,b,h,param,sex.rat,z5,bias)
{
  with(as.list(c(param,sex.rat)),
       {
         omega<-0
         ht<-1-h
         if(bias){
           Phimin<-rep(0,length(Sm))
           for(i in 1:length(Sm)){
             state <- c(z1=0.5,z3=0.5,z5=0.5,h=0,omega=0)
             parameters <- c(F1=F1,Sf=Sf,Sm=Sm[i],b=b,K=0.1,phi=0)
             times <- seq(0,30)
             out <- ode(y=state,times=times,func=Evolution,parms=parameters)
             z1<-out[length(out[,2]),2]
             z3<-out[length(out[,3]),3]
             if(z5>0){z5<-z5}
             else{z5<-z3}
             z1t<-1-z1
             z3t<-1-z3
             z5t<-1-z5
             F3<-F1
             omegat<-1-omega
             H<-h*z1t*F1
             pi<-1-exp(-H)
             pit<-1-pi
             F5<-F3*(1+b*H)
             B<-Sf*b
             u3<-ht*z1t*F1
             u5<-Sf
             G<-Sf*F3*pit+(1/2)*ht*z1t*F1*F3
             L<-Sf*F5*pi*z5*(omegat)
             Om<-(Sm[i]*z1*F1)/(z3*G+Sm[i]*z1*F1+L)
             Im<-Sf*F5*pi*((z3-z5)/z3t)/(z3*G+Sm[i]*z1*F1+L)
             Phimin[i]<-(z3t/z3)*(1-Om+Im)/(3*z5t+z5*((z3t/z3)*(1-Om+Im)))
           }
           return(Phimin)
         }
         else{
           z1t<-1-z1
           z3t<-1-z3
           z5t<-1-z5
           F3<-F1
           omegat<-1-omega
           H<-h*z1t*F1
           pi<-1-exp(-H)
           pit<-1-pi
           F5<-F3*(1+b*H)
           B<-Sf*b
           u3<-ht*z1t*F1
           u5<-Sf
           
           G<-Sf*F3*pit+(1/2)*ht*z1t*F1*F3
           
           L<-Sf*F5*pi*z5*(omegat)
           
           Om<-(Sm*z1*F1)/(z3*G+Sm*z1*F1+L)
           
           Im<-Sf*F5*pi*((z3-z5)/z3t)/(z3*G+Sm*z1*F1+L)
           
           Phim<-(z3t/z3)*(1-Om+Im)/(3*z5t+z5*((z3t/z3)*(1-Om+Im)))
           
           return(Phim)
         }
       })
}


# Coevolution function ---------------------------------------------------------

Evolution <- function(t,state,parameters)
{
  with(as.list(c(state,parameters)),
       {
         ## test
         
         #   h<-0.5
         #   omega<-0.5
         #   z1<-0.6
         #   z3<-0.5
         #   z5<-0.7
         #   F1<-5
         #   Sm<-0.5
         #   Sf<-0.8
         #   phi<-0.8
         #   b<-1.1
         #   K<-0.001
         
         
         # Expressions
         
         ht<-1-h
         z1t<-1-z1
         z3t<-1-z3
         z5t<-1-z5
         F3<-F1
         omegat<-1-omega
         H<-h*z1t*F1
         pi<-1-exp(-H)
         pit<-1-pi
         F5<-F3*(1+b*H)
         B<-Sf*b
         u3<-ht*z1t*F1
         u5<-Sf
         
         ZA<-pi*z5*(omegat+z5*omega*phi)+pit*z3
         
         ZAt<-pi*z5t*(1+z5*omega*phi)+pit*z3t
         
         V<-(1/2)*ht*z3t*F3
         
         G<-Sf*F5*pit+(1/2)*ht*z1t*F1*F3
         
         L<-Sf*F5*pi*z5*(omegat+z5*omega*phi)
         
         Om<-(Sm*z1*F1)/(z3*G+Sm*z1*F1+L)
         
         Im<-Sf*F5*pi*((1+z5*phi*omega)*((z3-z5)/z3t)+z5*omega)/(z3*G+Sm*z1*F1+L)
         
         (vf1<-V*(2-Om+Im))
         (vm1<-V*(z1t/z1)+(1/2)*Sm*(z3t/z3)*(1-Om+Im))
         (vm2<-(1/2)*(z3t/z3)*(1-Om+Im))
         
         ## Selection differentials
         
         
         (z1diff<-K*F1*(V*(-1+(z1t/z1)+(1/2)*Om-(1/2)*Im)+(1/2)*Sm*(z3t/z3)*(1-Om+Im)
                        -(1/2)*B*h*F3*(ZAt+(z3t/z3)*ZA*(1-Om+Im))))
         if (z1<0) {z1diff <- 0}
         if (z1>1)  {z1diff <- 0}
         
         
         (z3diff<-K*(1/2)*((u3*F3+u5*F5*pit)/(u3+u5))*(-1+(z3t/z3)*(1-Om+Im)))
         if (z3<0) {z3diff <- 0}
         if (z3>1) {z3diff <- 0}
         
         
         (z5diff<-K*((1/2)*Sf*F5*pi/(u3+u5))*(-1+omega*phi*(1-2*z5)+(z3t/z3)*(1-Om+Im)*(omegat+2*omega*phi*z5))  )
         if (z5<0) {z5diff <- 0}
         if (z5>1)  {z5diff <- 0}
         if (t<=1000) {z5diff<-z3diff}
         
         
         (Hdiff<-K*(1/2)*F5*(-z3t*(2-Om+Im)+(1/2)*B*(3*z5t*(1+z5*omega*phi)
                                                     +z5*(omegat+z5*omega*phi)*(z3t/z3)*(1-Om+Im))))
         if (h<0.0001) {Hdiff <- 0}
         if (h>1) {Hdiff <- 0}
         
         
         #(omediff<-K*(1/4)*Sf*F5*h*B*z5*(phi*(3*z5t+z5*(z3t/z3)*(1-Om+Im))-z5*(z3t/z3)*(1-Om+Im)))/(u3+u5)
         omediff<-K*(1/4)*Sf*F5*h*B*z5*(phi*(3*z5t+z5*(z3t/z3)*(1-Om+Im))-(z3t/z3)*(1-Om+Im))/(u3+u5)
         if (omega<0.0001)  {omediff <- 0 }
         if (omega>1)       {omediff <- 0 }
         
         
         list(c(z1diff,z3diff,z5diff,Hdiff,omediff))
       }) 
} # Closes Evolution function
ExponentialGrowth<-function(x,lambda){
   return (lambda*x)
}

Euler<-function(x, D, FUN,lambda){
   return(x+FUN(x,lambda)*D)
}

ExponentialSolution<-function(N0,t,lambda){
   return (N0*exp(t*lambda))
}

N0<-0.25
lambda<-0.3
plot(c(0,0),col="0",xlim=c(0,10),ylim=c(0,(N0*exp(lambda*10)+1)))
for (Deltat in c(0.05,0.1,0.2)){
   PointsToEstimate<-seq(0,10,by=Deltat)
   Iterations<-length(PointsToEstimate)-1
   Approx<-rep(0,Iterations+1)
   RealSol<-rep(0,Iterations+1)
   Approx[1]<-RealSol[1]<-N0
   for (i in 1:Iterations){
      RealSol[i+1]<-ExponentialSolution(N0,Deltat*i,lambda)		
      Approx[i+1]<-Euler(Approx[i],Deltat,ExponentialGrowth,lambda)		
   }
   if (Deltat==0.05) mycol="red"
   if (Deltat==0.1) mycol="blue"
   if (Deltat==0.2) mycol="green"
   points(RealSol~PointsToEstimate,col="black",type="l")
   points(Approx~PointsToEstimate,col=mycol,type="l")
}

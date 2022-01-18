ExponentialGrowth<-function(x,lambda){
   # FUN
   return (lambda*x)
}

# approximation
Euler<-function(x, D, FUN,lambda){
   # x = previous time step, D = step size, fun = approximation, 
   # lambda = input for FUN
   # recursive calls?
   return(x+FUN(x,lambda)*D)
}

# exact solution
ExponentialSolution<-function(N0,t,lambda){
   # N0 = initial pop, t = current time, lambda = growth rate
   return (N0*exp(t*lambda))
}

N0<-0.25
lambda<-0.3

plot(c(0,0),col="0",
     xlab="time", ylab="population size",
     xlim=c(0,10),ylim=c(0,(N0*exp(lambda*10)+1)))
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

plot(c(0,0),col="0",
     xlab="step size delta t", ylab="error at time t=10",
     xlim=c(0,0.25),ylim=c(0,0.6))
errors <- c(0)
for (Deltat in c(0.05,0.1,0.2)){
   PointsToEstimate<-seq(0,10,by=Deltat)
   Iterations<-length(PointsToEstimate)-1
   Approx<-rep(0,Iterations+1)
   RealSol<-rep(0,Iterations+1)
   Approx[1]<-RealSol[1]<-N0
   for (i in 1:Iterations){
      Approx[i+1]<-Euler(Approx[i],Deltat,ExponentialGrowth,lambda)		
   }
   errors <- append(errors,5.0213842-Approx[length(Approx)])
}
points(errors ~ c(0,0.05, 0.1, 0.2), col = "black", type = "l")





#### Repeat for logistic growth
LogisticGrowth<-function(x,lambda){
   # FUN
   return (lambda*x*(1-x))
}

# approximation
Euler<-function(x, D, FUN,lambda){
   # x = previous time step, D = step size, fun = approximation, 
   # lambda = input for FUN
   # recursive calls?
   return(x+FUN(x,lambda)*D)
}

# exact solution
LogisticSolution<-function(N0,t,lambda){
   # N0 = initial pop, t = current time, lambda = growth rate
   return ((N0*exp(t*lambda))/(1+N0*(exp(t*lambda)-1)))
}

N0<-0.25
lambda<-0.3

plot(c(0,0),col="0",
     xlab="time", ylab="population size",
     xlim=c(0,10),ylim=c(0,(N0*exp(10*lambda))/(1+N0*exp(10*lambda)-1)))
for (Deltat in c(0.05,0.1,0.2)){
   PointsToEstimate<-seq(0,10,by=Deltat)
   Iterations<-length(PointsToEstimate)-1
   Approx<-rep(0,Iterations+1)
   RealSol<-rep(0,Iterations+1)
   Approx[1]<-RealSol[1]<-N0
   for (i in 1:Iterations){
      RealSol[i+1]<-LogisticSolution(N0,Deltat*i,lambda)		
      Approx[i+1]<-Euler(Approx[i],Deltat,LogisticGrowth,lambda)		
   }
   print(Deltat)
   print(RealSol)
   print(Approx)
   if (Deltat==0.05) mycol="red"
   if (Deltat==0.1) mycol="blue"
   if (Deltat==0.2) mycol="green"
   points(RealSol~PointsToEstimate,col="black",type="l")
   points(Approx~PointsToEstimate,col=mycol,type="l")
}


# plot the error terms
plot(c(0,0),col="0",
     xlab="step size delta t", ylab="error at time t=10",
     xlim=c(0,0.25),ylim=c(-0.005,0.005))
errors <- c(0)
for (Deltat in c(0.05,0.1,0.2)){
   PointsToEstimate<-seq(0,10,by=Deltat)
   Iterations<-length(PointsToEstimate)-1
   Approx<-rep(0,Iterations+1)
   RealSol<-rep(0,Iterations+1)
   Approx[1]<-RealSol[1]<-N0
   for (i in 1:Iterations){
      Approx[i+1]<-Euler(Approx[i],Deltat,LogisticGrowth,lambda)		
   }
   errors <- append(errors,(N0*exp(10*lambda))/(1+N0*(exp(10*lambda)-1))-Approx[length(Approx)])
}
points(errors ~ c(0,0.05, 0.1, 0.2), col = "black", type = "l")

print(errors)

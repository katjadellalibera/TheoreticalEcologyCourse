## ExpGrowth.R 
## Stefano Allesina sallesina@uchicago.edu allesinalab.uchicago.edu
## Code for "Theoretical Ecology"

require(deSolve) 
## package for integrating numerically ODEs see 
##http://cran.r-project.org/web/packages/deSolve/

### Exponential Growth ### 

## This function takes a step in time
ExponentialGrowth<-function(t, state, parameters) { 
  X <- state[1] ## the first element is the density of X at time t
  r <- parameters[1] ## the first parameter is the growth rate
  dX<-r*X 
  return(list(dX)) ## for some reason, you have to return a list
}

## This function runs the model and produces the trajectory
RunExponentialGrowth<-function(MaxTime=10,GrowthRate=0.2,InitialX=1.0){
  times <- seq(0, MaxTime, by = 0.01)
  parameters <- c(r=GrowthRate)
  state <- c(ExponentialX=InitialX)
  out <- ode(y = state, times = times, func = ExponentialGrowth, parms = parameters)
  print(plot(out[,2]~out[,1],type="l",xlab="time",ylab="Density of X"))
  return(out)
}

# run the solver for exponential growth 
out <- RunExponentialGrowth()



### Logistic Growth ###

## This function takes a step in time
LogisticGrowth<-function(t, state, parameters) { 
  X <- state[1] ## the first element is the density of X at time t
  r <- parameters[1] ## the first parameter is the growth rate
  k <- parameters[2] ## the second parameter is the carrying capacity
  dX<-r*X*(1-X/k) 
  return(list(dX)) ## for some reason, you have to return a list
}


## This function runs the model and produces the trajectory
RunLogisticGrowth<-function(MaxTime = 10, GrowthRate = 0.5,
                            CarryingCapacity = 10, InitialX = 1.0){
  times <- seq(0, MaxTime, by = 0.01)
  parameters <- c(r=GrowthRate,k=CarryingCapacity)
  state <- c(LogisticX=InitialX)
  out <- ode(y = state, times = times, func = LogisticGrowth, parms = parameters)
  print(plot(out[,2]~out[,1],type="l",xlab="time",ylab="Density of X"))
  return(out)
}

# run the solver for exponential growth 
logisticOut <- RunLogisticGrowth()



### Allee Growth ###

## This function takes a step in time
AlleeGrowth<-function(t, state, parameters) { 
  X <- state[1] ## the first element is the density of X at time t
  u <- parameters[1] ## the first parameter is the upper limit of the growth rate
  v <- parameters[2] ## the second parameter is critical point
  c <- parameters[3] ## the third parameter is the density dependent death rate
  dX<-((u*X/(v+X)-c*X)*X)
  return(list(dX)) ## for some reason, you have to return a list
}


## This function runs the model and produces the trajectory
RunAlleeGrowth<-function(MaxTime = 10, UpperGrowthRate = 2,
                          CriticalPoint = 5, DeathRate = 0.1, InitialX = 1.0){
  times <- seq(0, MaxTime, by = 0.01)
  parameters <- c(u=UpperGrowthRate,v=CriticalPoint, c= DeathRate)
  state <- c(AlleeX=InitialX)
  out <- ode(y = state, times = times, func = AlleeGrowth, parms = parameters)
  print(plot(out[,2]~out[,1],type="l",xlab="time",ylab="Density of X"))
  return(out)
}

# run the solver for exponential growth 
AlleeOut <- RunAlleeGrowth()



### Levin Growth ###

## This function takes a step in time
LevinGrowth<-function(t, state, parameters) { 
  X <- state[1] ## the first element is the density of X at time t
  c <- parameters[1] ## the first parameter is the colonization rate
  e <- parameters[2] ## the second parameter is the extinction
  dX<-(c*X*(1-X)-e*X)
  return(list(dX)) ## for some reason, you have to return a list
}


## This function runs the model and produces the trajectory
RunLevinGrowth<-function(MaxTime = 10, ColonizationRate = 0.7,
                          ExtinctionRate = 0.1, InitialX = 0.1){
  times <- seq(0, MaxTime, by = 0.01)
  parameters <- c(c= ColonizationRate, e= ExtinctionRate)
  state <- c(LevinX=InitialX)
  out <- ode(y = state, times = times, func = LevinGrowth, parms = parameters)
  print(plot(out[,2]~out[,1],type="l",xlab="time",ylab="Density of X"))
  return(out)
}

# run the solver for exponential growth 
LevinOut <- RunLevinGrowth()


### Theta-Logistic Growth ###

## This function takes a step in time
ThetaLogisticGrowth<-function(t, state, parameters) { 
  X <- state[1] ## the first element is the density of X at time t
  r <- parameters[1] ## the first parameter is the growth rate
  k <- parameters[2] ## the second parameter is the carrying capacity
  theta <- parameters[3] ## the third parameter is the theta
  dX<-r*X*(1-(X/k)^theta)
  return(list(dX)) ## for some reason, you have to return a list
}


## This function runs the model and produces the trajectory
RunThetaLogisticGrowth<-function(MaxTime = 10, GrowthRate = 1,
                         CarryingCapacity = 10, Theta = 0.5,
                         InitialX = 1){
  times <- seq(0, MaxTime, by = 0.01)
  parameters <- c(r= GrowthRate, k = CarryingCapacity, theta = Theta)
  state <- c(ThetaX=InitialX)
  out <- ode(y = state, times = times, func = ThetaLogisticGrowth, parms = parameters)
  print(plot(out[,2]~out[,1],type="l",xlab="time",ylab="Density of X"))
  return(out)
}

# run the solver for exponential growth 
ThetaLogisticOut <- RunThetaLogisticGrowth()



### Combine into one plot
allOut <- merge(as.data.frame(out),as.data.frame(logisticOut), 
                by = "time")
allOut <- merge(allOut,as.data.frame(AlleeOut), by = "time")
allOut <- merge(allOut,as.data.frame(LevinOut), by = "time")
allOut <- merge(allOut,as.data.frame(ThetaLogisticOut), by = "time")


ggplot(allOut,  aes(x=time))  +
  geom_line(aes(y=ExponentialX, colour = "Exponential")) +
  geom_line(aes(y=LogisticX, colour = "Logistic")) +
  geom_line(aes(y=AlleeX, colour = "Allee")) +
  geom_line(aes(y=ThetaX, colour = "Theta-Logistic")) +
  scale_color_manual(name = "Model", 
                     values = c("Exponential" = "red", "Logistic" = "blue",
                                "Allee" = "orange", 
                                "Theta-Logistic" = "darkgreen")) +
  ylab("Population Size") +
  xlab("Time")


# plot Levin
ggplot(as.data.frame(LevinOut)) +
  aes(x=time, y=LevinX) +
  geom_line() +
  ylab("Fraction of Occupied patches") 

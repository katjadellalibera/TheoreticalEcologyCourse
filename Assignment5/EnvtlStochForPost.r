#Declaring parameter values
lambda = 2.0;
Rlzns = 10000;
StDev = 0.5;
TimeSteps = 10;
InitN = 1.0;

#Declaring that N and Stor are vectors
N <- numeric();
Stor <- numeric();

#When we plot N vs t, it is helpful to have t already holding the values of time
t = c(0:(TimeSteps));



#We are going to plot two graphs, the next line sets that up for us
par(mfrow=c(2,1));

Stor <- numeric(); # just putting this here so I don't forget to rerun it
for(j in 1:Rlzns){ #This is a for loop, so that we can plot many realizations
	N[1] = InitN  #1. Change this line so that N[1] gets the initial value of N
	for(i in 2:(TimeSteps+1)){ #This for loop iterates the model over the number of time steps. 
		N[i] = lambda*exp(rnorm(1,0,StDev))*N[i-1] #2. Change this line so that N[i] is the new value of the population size
	}
		
	Stor[j] = N[TimeSteps+1] #3. Change this line so that Store gets the LAST value of N

	if(j==1){ #For the first realization(visualization?), we use "plot"
			plot(t,log(N),type="l",xaxs="i",ylim=c(-3,12),xlim=c(0,10),xlab="Time (generations)",ylab="log Pop. Size");
	}else{ #For subsequent realizations, we use "lines", to add lines to the original plot
			lines(t,log(N));	
	}
}


#Next line is a function that plots a straight line, with intercept "a" and slope "b"
abline(a=0,b=log(lambda),col="RED",lwd=4); #4. Ca
	
#Plot a histogram of the population sizes in "Stor"
hist(log(Stor),freq=FALSE,breaks=200,right=TRUE,
     xlab="Final log population size",ylab="Frequency")
#Create a vector "x"...
x = seq(from=-20,to=20,by=0.1);
TheoryAvg = log(InitN) + TimeSteps*log(lambda) #5.  Insert the theoretical expectation expression here 
TheoryVar = TimeSteps*StDev^2#6. Insert the theoretical variance expression here
#Calculate the probability distribution function of a normal distribution using the theoretical expectation and variance
Dist = dnorm(x,mean=TheoryAvg,sd=TheoryVar^(0.5));
#Plot the distribution function
lines(x,Dist,lwd=1,col="RED");



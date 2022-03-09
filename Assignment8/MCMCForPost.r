#A function for plotting error bars
ErrBars<-function(x,y,SE){
	NumPts = length(x);
	for(Pts in 1:NumPts){
		xvals = c(x[Pts],x[Pts]);
		Upper = y[Pts]+SE[Pts];
		Lower = y[Pts]-SE[Pts];
		yvals = c(Upper,Lower);
		lines(xvals,yvals);
	}
}


#Reading in the data
Data = read.csv("VirusData6.csv");



######### Downhill Simplex to Generate Proposal Means

#The log-likelihood function with heterogeneity
logLHoodHtg<-function(par){

	betaBar = par[1];
	V  = par[2];

	timeT = 7.0;
	qoft = (1+betaBar*V*Data$Density*timeT)**(-1/V); #1. Type in the model that describes the probability of surviving to time T
	n0 = Data$Total; #n0 is the total number of insects in each treatment
	N = Data$Uninfected; #N is the number uninfected
	logLHood = log(choose(n0,N)*qoft**N*(1-qoft)**(n0-N)) ; #2. Type in the log-likelihood function for the pure death model
	return(-sum(logLHood));  #We are looking for the negative sum of the log likelihoods across experimental units, because optim is a minimizer

}



#Finding the best-fit parameter values
HtgOptOut = optim(par=c(0.1,0.1),fn=logLHoodHtg);
BestBeta = HtgOptOut$par[1]; 
BestVee = HtgOptOut$par[2];

Tmts = unique(Data$Density); #Figure out what the treatments were
AvgTransmn = tapply(-log(Data$Uninfected/Data$Total),Data$Density,mean); #Use tapply to calculate the mean transmission rateat each treatment
SDTransmn = tapply(-log(Data$Uninfected/Data$Total),Data$Density,sd); #Use tapply to calculate the standard deviation of the transmission rate

par(mfrow=c(1,1));
par(mai=c(1,1,1,1));
plot(Tmts,AvgTransmn,ylim=c(0,2),xlim=c(0,max(Tmts)),xlab="Cadavers per branch",ylab="Transmission"); #Plot the model prediction against the data


NumReps = nrow(Data)/length(Tmts); #Count how many replicates there were, for the standard-error calculation
ErrBars(x=Tmts,y=AvgTransmn,SDTransmn/sqrt(NumReps-1)); #Plot the error bars - standard error is sd/sqrt(n-1)


Cadavers = seq(from=0,to=4,by=0.1);
time = 7;
BestHtgModel = (1/BestVee)*log(1 + BestBeta*time*BestVee*Cadavers);
lines(Cadavers,BestHtgModel,col="RED",lwd=2.5);




############## MCMC: running the chains #####################

maxitn = 2e5; maxchain = 5;
LH = c(); MaxPars = 2;

ThinStop = 1000;
ThinStor = maxitn/ThinStop;
StorPars = array(dim=c(maxchain,ThinStor,MaxPars));

PropMean = c(); PropSD = c(); Prior = c()
#2. Fix the proposals.  Do the means on a log scale.
PropMean[1] = log(BestBeta); 
PropMean[2] = log(BestVee); 
PropSD[2] = 2.0;  
PropSD[1] = 2.0; 


Prior[1] = 1000.0;
Prior[2] = 1000.0;
ScaleSD = 1.2;


for(chain in 1:maxchain){
	CurrentPars = c(); #the current values of the parameters
	NewPars = c(); #the proposed values of the parameters
	for(i in 1:MaxPars){ #starting values for the parameters
		CurrentPars[i] = exp(rnorm(1,mean=PropMean[i],sd = ScaleSD*PropSD[i]));
	}
	#Calculating the posterior of the initial parameter set
	OldLH = -1.0*logLHoodHtg(par=CurrentPars); #-1 because logLHoodHtg returns -log lhood, and we want log lhood
	OldPost = OldLH;
	for(i2 in 1:MaxPars){
		OldPost = OldPost + log(dunif(x=(CurrentPars[i2]),max=Prior[i2]));
	}
	#parHtg = c(exp(BestBeta),exp(BestVee));
	#logLHoodHtg(parHtg);

	AcceptCount = 0;
	itnStor = 1;
	for(itn in 1:maxitn){
		if(itn==1){
			ScaleSD = 1.2;
		}else{
			ScaleSD = 1.0;
		}
		for(i in 1:MaxPars){
			for(i2 in 1:MaxPars){
				 NewPars[i2] = CurrentPars[i2];
			}
			#3.  Propose a new parameter value, and calculate its posterior.
			NewPars[i] = exp(rnorm(1, mean = PropMean[i], sd = PropSD[i]));
			NewLH = -1*logLHoodHtg(par=NewPars); 
			NewPost = NewLH;
			for(i2 in 1:MaxPars){
			 	NewPost = NewPost + log(dunif(x=(NewPars[i2]),max=Prior[i2]));
	      	}
			
			#Here is the adjustment for the proposal for the current parameter value
			OldPropAdj = log(dnorm(x=log(CurrentPars[i]), mean = PropMean[i],sd=ScaleSD*PropSD[i]));
			
			
			#4.  Now do the adjustment for the proposal, based on the new parameter value
			NewPropAdj = log(dnorm(x=log(NewPars[i]), mean = PropMean[i], sd=PropSD[i]))

			
			if(itn>1){
			 	#6. Write out the acceptance criterion, and decide whether to accept or not 
				#6a. First, calculate the acceptance criterion
				Criterion = exp(NewPost - NewPropAdj - (OldPost - OldPropAdj));
				#6b. Second, draw a U(0,1) random variate.
			 	temp = runif(1,0,1);
				#6c. If statement for whether or not to accept
			 	if(Criterion>temp){
					AcceptCount = AcceptCount + 1;
					OldPost = NewPost;
					OldLH = NewLH;
					#7. Acceptance means setting current parameters equal to proposed parameters
					for(i2 in 1:MaxPars){
				 		CurrentPars[i2] = NewPars[i2];
					}#for i2
			 	} #if Criterion
			} #itn>1
		} #for i in 1:MaxPars

	    #We only want to keep some of the iterations
	    if((itn %% ThinStop)==0){
		cat("chain:",chain,"itn:",itn,"itnStor:",itnStor);
		for(i2 in 1:MaxPars){
			StorPars[chain,itnStor,i2] = CurrentPars[i2];
			cat(" i2:",i2,"StorPars:",StorPars[chain,itnStor,i2]);
		}
		itnStor = itnStor + 1;
		cat("\n"); #line = readline();
		
	    } #Thin

					
	} #itn

} #chain



########### MCMC diagnostics: did the algorithm converge? #############
require(coda);
AllPars = c();
itnStor = itnStor - 1;
HalfWay = round(0.5*itnStor);
x1 = mcmc(StorPars[1,HalfWay:itnStor,]);
x2 = mcmc(StorPars[2,HalfWay:itnStor,]);
x3 = mcmc(StorPars[3,HalfWay:itnStor,]);
x4 = mcmc(StorPars[4,HalfWay:itnStor,]);
x5 = mcmc(StorPars[5,HalfWay:itnStor,]);
xAll = mcmc.list(x1,x2,x3,x4,x5);

#Trace plots
par(ask="TRUE");
par(mai=c(0.5,0.5,0.5,0.5));
plot(xAll)
gelmanOut = gelman.diag(xAll)
print(gelmanOut);

#Autocorrelations
par(mfrow=c(5,2));
par(ask="TRUE");
par(mai=c(0.225,0.25,0.225,0.25));
#par(mai=c(0.5,0.5,0.5,0.5));
title(main="")
acf(x1[,1],xlab="",ylab=""); abline(h=-0.1,col="RED"); abline(h=0.1,col="RED");
acf(x1[,2],xlab="",ylab=""); abline(h=-0.1,col="RED"); abline(h=0.1,col="RED");
acf(x2[,1],xlab="",ylab=""); abline(h=-0.1,col="RED"); abline(h=0.1,col="RED");
acf(x2[,2],xlab="",ylab=""); abline(h=-0.1,col="RED"); abline(h=0.1,col="RED");
acf(x3[,1],xlab="",ylab=""); abline(h=-0.1,col="RED"); abline(h=0.1,col="RED");
acf(x3[,2],xlab="",ylab=""); abline(h=-0.1,col="RED"); abline(h=0.1,col="RED");
acf(x4[,1],xlab="",ylab=""); abline(h=-0.1,col="RED"); abline(h=0.1,col="RED");
acf(x4[,2],xlab="",ylab=""); abline(h=-0.1,col="RED"); abline(h=0.1,col="RED");
acf(x5[,1],xlab="",ylab=""); abline(h=-0.1,col="RED"); abline(h=0.1,col="RED");
acf(x5[,2],xlab="",ylab=""); abline(h=-0.1,col="RED"); abline(h=0.1,col="RED");

	

############### Decision-making: should we spray the virus? #################
AllData = c();
for(chain in 1:maxchain){
	AllData = rbind(AllData,StorPars[chain,1:ThinStor,]);
}


Burnout<-function(nu,mu,V,S0,I0,i){ #This function is a straight line

	return (1 - i - (1 + nu*V/mu*(i*S0+I0))^(-V));

}

S0 = 0.1; I0 = 0.01;
UpperPctl = round(0.25*ThinStor);
LowerPctl = round(0.75*ThinStor);
S0Range = seq(from=0.0,to=10.0,by=0.1); #the range of values over which we will calculate the root
S0Dim = length(S0Range);
ParamDim = nrow(AllData);
#FractInf = array(dim=c(S0Dim,ParamDim));
FractInf = c();
FIPlotUpper = c(); FIPlotLower = c(); FIMedian = c();
i = 1;
for(S0 in S0Range){
	for(j in 1:ThinStor){
	
		z = uniroot(f=Burnout,interval=c(0.0,1.0),nu=AllData[j,1],mu=0.41,V=AllData[j,2],S0=S0,I0=I0);
		FractInf[j] = z$root;
		#cat("S0:",S0,"j:",j,"FractInf:",FractInf[j],"\n");

		if(1){
			Ordering = order(FractInf,decreasing=TRUE); #figuring out the order of the parameters
	 		OrderedFractInf = cbind(FractInf)[Ordering,] #Carrying out the ordering...
			FIPlotUpper[i] = OrderedFractInf[UpperPctl];
			FIPlotLower[i] = OrderedFractInf[LowerPctl];
			FIMedian[i] = median(FractInf);
		}

	}
	#line = readline();
	i = i + 1;
}

par(ask="FALSE");



k1 = 0.0;
k2 = 100.0;
#k2 = 0.05;

par(mfrow=c(1,1));

MedianCost = k1 + k2*(1-FIMedian)*S0Range;
LowerCost = k1 + k2*(1-FIPlotUpper)*S0Range;
UpperCost = k1 + k2*(1-FIPlotLower)*S0Range;
#plot(S0Range,(1-FIMedian)*S0Range*k1*exp(k2*Defoln),type="l"); #exponential
par(mfrow=c(1,1));
Upperylim = max(UpperCost);
Lowerylim = min(LowerCost);
par(mai=c(1,1,1,1));
plot(S0Range,MedianCost,type="l",ylim=c(Lowerylim,Upperylim),ylab="Cost",xlab="Insect Density");  #linear
lines(S0Range,UpperCost,type="l",lwd=2.5);
lines(S0Range,LowerCost,type="l",lwd=2.5);










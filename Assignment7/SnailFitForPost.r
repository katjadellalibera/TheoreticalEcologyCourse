
#reading in the Data
Data = read.csv("GregsVer4.csv");


########################PLOTTER()###########################################
#Plotter() is a function that plots the model predictions against the data.
#It is a bit complicated, so I wrote it for you.
#We need it because whenever you use a nonlinear fitting routine,
#you should always visually compare the model to the data, if you can.
############################################################################
Plotter<-function(params){
	MaxPlots = length(unique(Data$Plot));
	x = numeric(); y = numeric();
	m = params[1]; I = params[2]; alpha = params[3];
	for(plt in 1:MaxPlots){

		Data2 = Data[Data$Plot==plt,];
		
		PredNoPred = "";
		if(Data2$Tmt[1]=="Pred"){
			mhat = m + alpha;
			if(plt==2){
				PredNoPred = "Predator";
			}
		}else{
			mhat = m;
			if(plt==8){
				PredNoPred = "No Predator";
			}

		}
		plot(Data2$Time,Data2$N,ylim=c(0,1200),xlab="Time",ylab="Number",main=PredNoPred,pch=19);
		Theory = Model(N=Data2,pars=c(mhat,I));
		#lines(Data3$Time,Theory);
		Times = unique(Data2$Time);
		tnum = length(Times) - 1; #subtract 1 becuz no theory prediction at last time point
		for(t in 1:tnum){
			y[1] = Data2$N[t];
			y[2] = Theory[t];
			x[1] = Times[t];
			x[2] = Times[t+1];
			lines(x,y);
			points(x,y,cex=2.0);
		}
	}
	return(1);
}

#####################End of Plotter()##########################################



########## Part 1. Fix up the model function.
Model<-function(N,pars){
	m = pars[1]; I = pars[2];
	Model = numeric();
	
	maxt = length(N$Time) - 1;
	for(t in 1:maxt){
		TimeDiff = N$Time[t+1]-N$Time[t]; #1.1 This is the difference between one sampling time and the next
		Model[t] = N$N[t]*exp(-m*TimeDiff)+I*(1-exp(-m*TimeDiff))/m; #1.2 The model goes here
	}

	
	return(Model);
}

########## Task 2. Fix up the code for the likelihood function without predation.
LHoodNoPred<-function(pars){

	MaxPlots = length(unique(Data$Plot));
	SSE = 0.0; 
	for(plt in 1:MaxPlots){

		Data2 = Data[Data$Plot==plt,];
		Theory = Model(N=Data2,pars=pars);
		Data3 = Data2[2:nrow(Data2),];
		SSE = SSE + sum((Data3$N-Theory)**2);  #2.1 Here you update the sum of squared errors between the data and the model

	}
	return(SSE);
}

############## Task 3. Fix up the code for the likelihood function with predation 
LHoodPred<-function(params){

	MaxPlots = length(unique(Data$Plot));
	SSE = 0.0; 
	m = params[1]; I = params[2]; alpha = params[3];
	for(plt in 1:MaxPlots){

		Data2 = Data[Data$Plot==plt,];
		if(Data2$Tmt[1]=="Pred"){
			mhat = m + alpha; #3.1 This is the death rate in the predation plots, which includes alpha
		}else{
			mhat = m; #3.2 This is the death rate in the predator-exclusion plots, which does not include alpha
		}
		pars = c(mhat,I); #Model parameters, to pass to the Model() function
		Theory = Model(N=Data2, pars = pars); #3.3 Calculate the model prediction here
		Data3 = Data2[2:nrow(Data2),];
		SSE = SSE + sum((Data3$N-Theory)**2); #3.4 Update the sum of squared errors

	}
	return(SSE);
}


################# Task 4. Write an AIC function
AICCalc<-function(n,K,LHood){
	return(n*log(LHood/n)+2*K*(n/(n-K-1))); #4.1 Type in the AIC equation here

}

################# Task 5.  Using optim() and finishing the calculation
#Now we calculate the AICc, which is n*log(SSE/n) + 2*K*n/(n-K-1)
#First we calculate n, then we calculate the likelihoods, and we use the 2 to calc. the AICc
#The initial population sizes don't count as samples, so n is the number of samples minus the inits
InitPops = length(Data$Time[Data$Time==0]);
n = nrow(Data) - InitPops; 

#Calculate the sum of squares for the no-predation model
Output = optim(par = c(1e-4,1.3),fn=LHoodNoPred);
#Calculate the AIC for the no-predation model
K1 = length(Output$par);
AIC1 = AICCalc(n=n,K=K1,LHood=Output$value) #AICCalc() is the function that calculates the AIC

#Calculate the sum of squares for the predation model
Output2 = optim(par = c(1e-4, 1.3,1), fn=LHoodPred); #5.1 Use the calculation for the no-predation model as a template here
#Calculate the AIC for the predation model
K2 = length(Output2$par);
AIC2 = AICCalc(n=n, K=K2, LHood=Output2$value); #5.2 Again use the calculation for the no-predation model as a template here


#Print out the AIC scores for each model, along with the AIC difference
cat("AIC1:",AIC1," AIC2:",AIC2," Delta AIC:",AIC1-AIC2,"\n");

#plot the data and the model prediction
par(mfrow=c(2,3));
par(ask=1);
Plotter(params=Output2$par);




## Thoughts: The AIC scores, (regardless of your choice of alpha) reflect some uncertainty 
# no matter what. We can't know the actual likelihood function.



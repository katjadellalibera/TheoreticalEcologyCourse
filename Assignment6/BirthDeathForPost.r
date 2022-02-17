#Parameters
MaxRlzns = 50; #NOTICE THAT THIS IS A SMALL NUMBER OF REALIZATIONS
NumExtinct = 0;
InitN = 2;
lambda = 6.0;
mu = 4.0;
TimePoints = c(0.1,0.5,1,2,3); #lambda = 6.0, mu = 3.6
NumTimePoints = length(TimePoints);
NStor <- array(0,c(NumTimePoints,MaxRlzns));
par(mfrow=c(2,1));
par(mai = c(0.85,0.8,0.3,0.25));#bottom,left,top,right


for(Rlzn in 1:MaxRlzns){

	N = numeric();
	t = numeric();

	N[1] = InitN;
	t[1] = 0.0;
	i = 1;
	Point = 1;
	while((N[i]>0)&&(t[i]<=TimePoints[NumTimePoints])){

		#################################
		# I. The first chunk of code that you must fix starts here: the simulation
		#################################
		AvgTimeToNextEvent = (lambda+mu)*N[i]; #1. Fill in the average time to the next event
		NextTime = rexp(1,AvgTimeToNextEvent); 
		i = i + 1;
		t[i] = t[i-1]+ NextTime; #2. Update the time by adding NextTime to the previous value of t (look at the preceding line if you are confused)

		rand = runif(1); #3. Here you must draw a single, uniformly distributed random variate, between 0 and 1
		BirthRate = lambda/(lambda+mu); #4.  This is the probability of a birth
		if(rand<BirthRate){ #5. If rand is less than the probability of a birth...
			N[i] = N[i-1]+1; #6. ...then a birth occurs, otherwise...
		} else{
			N[i] = N[i-1]-1; #7. ...a death occurs
		}
		if(N[i]<1) NumExtinct = NumExtinct+1; #8.  If the population has gone extinct, increase the number of extinctions by 1
			
		stop = 0;
		while((Point<=NumTimePoints)&&(stop!=1)){ #This while statement saves the population size at predetermined time points
			if(t[i]>TimePoints[Point]){
				
				NStor[Point,Rlzn] = N[i-1];
				Point = Point + 1;
			} else { stop = 1;}

		}
		
		#i = i + 1;
		

	} # while loop

	if(Rlzn==1){
		plot(t,log(N),type="l",xlab="Time",ylab="log Population Size",ylim=c(0,9),xlim=c(0,TimePoints[NumTimePoints]));
	}else{
		if(Rlzn<=20){
			if(min(N)<1){
				lines(t,log(N),col="RED");
			}else{
				lines(t,log(N));
			}
		}

	}
	
} # rlzn loop

#############################
#II. The second piece of code that you must fix starts here: the theory calculations
#############################
TheoryProb = (mu/lambda)^InitN; #9. This should the theoretical probability of extinction
ObsvdProb = NumExtinct/MaxRlzns; #10. This should be the observed frequency of extinctions 
PctError = (TheoryProb - ObsvdProb)/ObsvdProb; 
TheoryCV = sqrt((lambda+mu)/(InitN*(lambda-mu))); #11. This should be the theoretical, long-term C.V. from lecture, for the lambda > mu case


#Graphics stuff
TheoryName = "Theor. Ext. Prob.:";
TheoryLgnd = paste(TheoryName,round(TheoryProb*100.0)/100.0);
ObsvdName = "Obsvd. Ext. Prob.:";
ObsvdLgnd = paste(ObsvdName,ObsvdProb);
text(1,7,TheoryLgnd);
text(1,6,ObsvdLgnd);

#Printing out the theoretical and observed probabilities of extinction, along with the percent error
cat("TheoryProb:",round(100*TheoryProb)/100," ObsvdProb:",ObsvdProb," PctError:",PctError,"\n");

#Calculating the observed C.V.
ObsvdCV = numeric();
for(Point in 1:NumTimePoints){
	ObsvdCV[Point] = sd(NStor[Point,])/mean(NStor[Point,]);

}


#The second graph plots the observed and theoretical C.V.
yHi = max(TheoryCV,ObsvdCV);
cat("TheoryCV:",TheoryCV," ObsvdCV:",ObsvdCV,"\n");
plot(TimePoints,ObsvdCV,type="o",ylim=c(0,yHi),xlab="Time", ylab = "C.V. of Pop. Size");
abline(h=TheoryCV,lwd=2);



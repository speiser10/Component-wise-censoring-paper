

#preventable interval censoring: dual censoring simulation study 
#jaime lynn speiser
#this is the code for the midpoint of the interval analysis

#load r libraries
library(simsurv)
library(survival)
library(icenReg)

#suppress warnings
options(warn = -1)

#parameters to input
###n: sample size
###beta_right: treatment beta coefficient for death (right censored outcome)
###beta_int: treatment beta coefficient for dementia/persistent mobility (interval censored outcome)
###n_ints: number of intervals
###t: total time of followup (years)

#right censored: simulated from simsurv package using weibull proportional hazards model
#interval censored: simulated based on function from icenReg package (edited) using weibull proportional hazards

#functions 
#edited from simIC_weib function from icenReg package
simsurv_Int<-function (n = 100, b1=-0.5, trt=trt, model = "ph", shape = 2, 
    scale = 2, inspections = 2, inspectLength = 2.5, rndDigits = NULL, 
    prob_cen = 1) {
    rawQ <- runif(n)
    nu <- exp(trt * b1)
    if (model == "ph") 
        adjFun <- function(x, nu) {
            1 - x^(1/nu)
        }
    else if (model == "po") 
        adjFun <- function(x, nu) {
            1 - x * (1/nu)/(x * 1/nu - x + 1)
        }
    adjQ <- adjFun(rawQ, nu)
    trueTimes <- qweibull(adjQ, shape = shape, scale = scale)
    obsTimes <- runif(n = n, max = inspectLength)
    if (!is.null(rndDigits)) 
        obsTimes <- round(obsTimes, rndDigits)
    l <- rep(0, n)
    u <- rep(0, n)
    caught <- trueTimes < obsTimes
    u[caught] <- obsTimes[caught]
    l[!caught] <- obsTimes[!caught]
    if (inspections > 1) {
        for (i in 2:inspections) {
            oldObsTimes <- obsTimes
            obsTimes <- oldObsTimes + runif(n, max = inspectLength)
            if (!is.null(rndDigits)) 
                obsTimes <- round(obsTimes, rndDigits)
            caught <- trueTimes >= oldObsTimes & trueTimes < 
                obsTimes
            needsCatch <- trueTimes > obsTimes
            u[caught] <- obsTimes[caught]
            l[needsCatch] <- obsTimes[needsCatch]
        }
    }
    else {
        needsCatch <- !caught
    }
    u[needsCatch] <- Inf
    if (sum(l > u) > 0) 
        stop("warning: l > u! Bug in code")
    isCensored <- rbinom(n = n, size = 1, prob = prob_cen) == 
        1
    l[!isCensored] <- trueTimes[!isCensored]
    u[!isCensored] <- trueTimes[!isCensored]
    if (sum(l == Inf) > 0) {
        allTimes <- c(l, u)
        allFiniteTimes <- allTimes[allTimes < Inf]
        maxFiniteTime <- max(allFiniteTimes)
        l[l == Inf] <- maxFiniteTime
    }
    return(data.frame(l = l, u = u))
}





#initialize parameters
n<-beta_right<-t<-n_ints<-beta_int<-NA

#simulation for using upper value of the interval
simRun2<-function(n,beta_right,beta_int,lambda_right,gammas_right,shape_int,scale_int){
	#create trt variable and time to death data (right censored)
	trt = rbinom(n, 1, 0.5)
	cov <- data.frame(id = 1:n, trt)
	#generate followup time, t, from random uniform between 3 and 5
	t<-runif(n,3,5)
	dat1 <- simsurv(lambdas = lambda_right, gammas = gammas_right, betas = c(trt = beta_right),x = cov, maxt = t)
	dat1 <- merge(cov,dat1)

	#create time to dementia/mobility data (interval censored)
	dat2 <- simsurv_Int(n=n, b1=beta_int, trt=trt, model="ph",shape=shape_int, scale=scale_int, inspections=t, inspectLength=1, prob_cen=1)

	#calculate which time to event comes first/make that the outcome and make censoring indicator
	dat<-data.frame(cbind(dat1,dat2))
	dat$status_int<-ifelse(dat$u<10000,1,0)
	dat$cens<-ifelse((dat$status_int+dat$status)>0,1,0)
	#use midpoint of the interval simulated
	for(i in 1:n){
		dat$y[i]<-min(dat$eventtime[i],(dat$u[i]+dat$l[i])/2)
	}

	#run the Cox regression model
	dat$SurvObj<-with(dat,Surv(y,cens==1))
	coxmodel<-coxph(SurvObj~trt,data=dat)
	#print(coxmodel)

	#save beta_trt estimate, SD, CI, events info
	beta_trt_est<-coxmodel$coefficients[1]
	beta_trt_sd<-sqrt(coxmodel$var)
	beta_trt_cilow<-beta_trt_est+ qnorm(.025) * beta_trt_sd
	beta_trt_ciup<-beta_trt_est+ qnorm(.975) * beta_trt_sd
	beta_pval<-summary(coxmodel)$waldtest[3]
	eventrate<-coxmodel$nevent/n
	eventrate_right<-mean(dat$status)
	eventrate_int<-mean(dat$status_int)

	#return bias, 95% coverage probability, power, event rates
	return(c(bias_int = beta_trt_est - (beta_int), bias_right = beta_trt_est - (beta_right), 
		coverage_int = ((beta_int > beta_trt_cilow) && (beta_int < beta_trt_ciup)),
		coverage_right = ((beta_right > beta_trt_cilow) && (beta_right < beta_trt_ciup)),
		power=(beta_pval<0.05),event_rate=eventrate, event_rate_right=eventrate_right,
		event_rate_int=eventrate_int))	
}


######################################################################################

set.seed(654)
nsims<-1000
sim1<-rowMeans(replicate(nsims,simRun2(n=20000,beta_right=-0.3,beta_int=-0.3,
	lambda_right=0.03,gammas_right=0.6,shape_int=2,scale_int=7)))
sim2<-rowMeans(replicate(nsims,simRun2(n=20000,beta_right=-0.3,beta_int=-0.15,
	lambda_right=0.03,gammas_right=0.6,shape_int=2,scale_int=7)))
sim3<-rowMeans(replicate(nsims,simRun2(n=20000,beta_right=-0.3,beta_int=0,
	lambda_right=0.03,gammas_right=0.6,shape_int=2,scale_int=7)))
sim4<-rowMeans(replicate(nsims,simRun2(n=20000,beta_right=-0.3,beta_int=0.15,
	lambda_right=0.03,gammas_right=0.6,shape_int=2,scale_int=7)))
sim5<-rowMeans(replicate(nsims,simRun2(n=20000,beta_right=-0.3,beta_int=0.3,
	lambda_right=0.03,gammas_right=0.6,shape_int=2,scale_int=7)))

sim6<-rowMeans(replicate(nsims,simRun2(n=20000,beta_right=-0.15,beta_int=-0.3,
	lambda_right=0.03,gammas_right=0.6,shape_int=2,scale_int=7)))
sim7<-rowMeans(replicate(nsims,simRun2(n=20000,beta_right=-0.15,beta_int=-0.15,
	lambda_right=0.03,gammas_right=0.6,shape_int=2,scale_int=7)))
sim8<-rowMeans(replicate(nsims,simRun2(n=20000,beta_right=-0.15,beta_int=0,
	lambda_right=0.03,gammas_right=0.6,shape_int=2,scale_int=7)))
sim9<-rowMeans(replicate(nsims,simRun2(n=20000,beta_right=-0.15,beta_int=0.15,
	lambda_right=0.03,gammas_right=0.6,shape_int=2,scale_int=7)))
sim10<-rowMeans(replicate(nsims,simRun2(n=20000,beta_right=-0.15,beta_int=0.3,
	lambda_right=0.03,gammas_right=0.6,shape_int=2,scale_int=7)))

sim11<-rowMeans(replicate(nsims,simRun2(n=20000,beta_right=0,beta_int=-0.3,
	lambda_right=0.03,gammas_right=0.6,shape_int=2,scale_int=7)))
sim12<-rowMeans(replicate(nsims,simRun2(n=20000,beta_right=0,beta_int=-0.15,
	lambda_right=0.03,gammas_right=0.6,shape_int=2,scale_int=7)))
sim13<-rowMeans(replicate(nsims,simRun2(n=20000,beta_right=0,beta_int=0,
	lambda_right=0.03,gammas_right=0.6,shape_int=2,scale_int=7)))
sim14<-rowMeans(replicate(nsims,simRun2(n=20000,beta_right=0,beta_int=0.15,
	lambda_right=0.03,gammas_right=0.6,shape_int=2,scale_int=7)))
sim15<-rowMeans(replicate(nsims,simRun2(n=20000,beta_right=0,beta_int=0.3,
	lambda_right=0.03,gammas_right=0.6,shape_int=2,scale_int=7)))

sim16<-rowMeans(replicate(nsims,simRun2(n=20000,beta_right=0.15,beta_int=-0.3,
	lambda_right=0.03,gammas_right=0.6,shape_int=2,scale_int=7)))
sim17<-rowMeans(replicate(nsims,simRun2(n=20000,beta_right=0.15,beta_int=-0.15,
	lambda_right=0.03,gammas_right=0.6,shape_int=2,scale_int=7)))
sim18<-rowMeans(replicate(nsims,simRun2(n=20000,beta_right=0.15,beta_int=0,
	lambda_right=0.03,gammas_right=0.6,shape_int=2,scale_int=7)))
sim19<-rowMeans(replicate(nsims,simRun2(n=20000,beta_right=0.15,beta_int=0.15,
	lambda_right=0.03,gammas_right=0.6,shape_int=2,scale_int=7)))
sim20<-rowMeans(replicate(nsims,simRun2(n=20000,beta_right=0.15,beta_int=0.3,
	lambda_right=0.03,gammas_right=0.6,shape_int=2,scale_int=7)))

sim21<-rowMeans(replicate(nsims,simRun2(n=20000,beta_right=0.3,beta_int=-0.3,
	lambda_right=0.03,gammas_right=0.6,shape_int=2,scale_int=7)))
sim22<-rowMeans(replicate(nsims,simRun2(n=20000,beta_right=0.3,beta_int=-0.15,
	lambda_right=0.03,gammas_right=0.6,shape_int=2,scale_int=7)))
sim23<-rowMeans(replicate(nsims,simRun2(n=20000,beta_right=0.3,beta_int=0,
	lambda_right=0.03,gammas_right=0.6,shape_int=2,scale_int=7)))
sim24<-rowMeans(replicate(nsims,simRun2(n=20000,beta_right=0.3,beta_int=0.15,
	lambda_right=0.03,gammas_right=0.6,shape_int=2,scale_int=7)))
sim25<-rowMeans(replicate(nsims,simRun2(n=20000,beta_right=0.3,beta_int=0.3,
	lambda_right=0.03,gammas_right=0.6,shape_int=2,scale_int=7)))

table2a<-rbind(sim1,sim2,sim3,sim4,sim5,sim6,sim7,sim8,sim9,sim10,sim11,sim12,sim13,sim14,sim15,
	sim16,sim17,sim18,sim19,sim20,sim21,sim22,sim23,sim24,sim25)
table2<-cbind(table2a[,1],table2a[,3],table2a[,2],table2a[,4],table2a[,5])
colnames(table2)<-c("BiasBint","CoverageBint","BiasBright","CoverageBright","Power")

write.csv(table2,"Simulation N=20000 12Pct Midpoint 5050 Results.csv")
write.csv(table2a,"Simulation N=20000 12Pct Midpoint 5050 Results Full.csv")


####################################################################################

set.seed(654)
nsims<-1000
sim1<-rowMeans(replicate(nsims,simRun2(n=20000,beta_right=-0.3,beta_int=-0.3,
	lambda_right=0.03,gammas_right=0.8,shape_int=2,scale_int=10)))
sim2<-rowMeans(replicate(nsims,simRun2(n=20000,beta_right=-0.3,beta_int=-0.15,
	lambda_right=0.03,gammas_right=0.8,shape_int=2,scale_int=10)))
sim3<-rowMeans(replicate(nsims,simRun2(n=20000,beta_right=-0.3,beta_int=0,
	lambda_right=0.03,gammas_right=0.8,shape_int=2,scale_int=10)))
sim4<-rowMeans(replicate(nsims,simRun2(n=20000,beta_right=-0.3,beta_int=0.15,
	lambda_right=0.03,gammas_right=0.8,shape_int=2,scale_int=10)))
sim5<-rowMeans(replicate(nsims,simRun2(n=20000,beta_right=-0.3,beta_int=0.3,
	lambda_right=0.03,gammas_right=0.8,shape_int=2,scale_int=10)))

sim6<-rowMeans(replicate(nsims,simRun2(n=20000,beta_right=-0.15,beta_int=-0.3,
	lambda_right=0.03,gammas_right=0.8,shape_int=2,scale_int=10)))
sim7<-rowMeans(replicate(nsims,simRun2(n=20000,beta_right=-0.15,beta_int=-0.15,
	lambda_right=0.03,gammas_right=0.8,shape_int=2,scale_int=10)))
sim8<-rowMeans(replicate(nsims,simRun2(n=20000,beta_right=-0.15,beta_int=0,
	lambda_right=0.03,gammas_right=0.8,shape_int=2,scale_int=10)))
sim9<-rowMeans(replicate(nsims,simRun2(n=20000,beta_right=-0.15,beta_int=0.15,
	lambda_right=0.03,gammas_right=0.8,shape_int=2,scale_int=10)))
sim10<-rowMeans(replicate(nsims,simRun2(n=20000,beta_right=-0.15,beta_int=0.3,
	lambda_right=0.03,gammas_right=0.8,shape_int=2,scale_int=10)))

sim11<-rowMeans(replicate(nsims,simRun2(n=20000,beta_right=0,beta_int=-0.3,
	lambda_right=0.03,gammas_right=0.8,shape_int=2,scale_int=10)))
sim12<-rowMeans(replicate(nsims,simRun2(n=20000,beta_right=0,beta_int=-0.15,
	lambda_right=0.03,gammas_right=0.8,shape_int=2,scale_int=10)))
sim13<-rowMeans(replicate(nsims,simRun2(n=20000,beta_right=0,beta_int=0,
	lambda_right=0.03,gammas_right=0.8,shape_int=2,scale_int=10)))
sim14<-rowMeans(replicate(nsims,simRun2(n=20000,beta_right=0,beta_int=0.15,
	lambda_right=0.03,gammas_right=0.8,shape_int=2,scale_int=10)))
sim15<-rowMeans(replicate(nsims,simRun2(n=20000,beta_right=0,beta_int=0.3,
	lambda_right=0.03,gammas_right=0.8,shape_int=2,scale_int=10)))

sim16<-rowMeans(replicate(nsims,simRun2(n=20000,beta_right=0.15,beta_int=-0.3,
	lambda_right=0.03,gammas_right=0.8,shape_int=2,scale_int=10)))
sim17<-rowMeans(replicate(nsims,simRun2(n=20000,beta_right=0.15,beta_int=-0.15,
	lambda_right=0.03,gammas_right=0.8,shape_int=2,scale_int=10)))
sim18<-rowMeans(replicate(nsims,simRun2(n=20000,beta_right=0.15,beta_int=0,
	lambda_right=0.03,gammas_right=0.8,shape_int=2,scale_int=10)))
sim19<-rowMeans(replicate(nsims,simRun2(n=20000,beta_right=0.15,beta_int=0.15,
	lambda_right=0.03,gammas_right=0.8,shape_int=2,scale_int=10)))
sim20<-rowMeans(replicate(nsims,simRun2(n=20000,beta_right=0.15,beta_int=0.3,
	lambda_right=0.03,gammas_right=0.8,shape_int=2,scale_int=10)))

sim21<-rowMeans(replicate(nsims,simRun2(n=20000,beta_right=0.3,beta_int=-0.3,
	lambda_right=0.03,gammas_right=0.8,shape_int=2,scale_int=10)))
sim22<-rowMeans(replicate(nsims,simRun2(n=20000,beta_right=0.3,beta_int=-0.15,
	lambda_right=0.03,gammas_right=0.8,shape_int=2,scale_int=10)))
sim23<-rowMeans(replicate(nsims,simRun2(n=20000,beta_right=0.3,beta_int=0,
	lambda_right=0.03,gammas_right=0.8,shape_int=2,scale_int=10)))
sim24<-rowMeans(replicate(nsims,simRun2(n=20000,beta_right=0.3,beta_int=0.15,
	lambda_right=0.03,gammas_right=0.8,shape_int=2,scale_int=10)))
sim25<-rowMeans(replicate(nsims,simRun2(n=20000,beta_right=0.3,beta_int=0.3,
	lambda_right=0.03,gammas_right=0.8,shape_int=2,scale_int=10)))

table2a<-rbind(sim1,sim2,sim3,sim4,sim5,sim6,sim7,sim8,sim9,sim10,sim11,sim12,sim13,sim14,sim15,
	sim16,sim17,sim18,sim19,sim20,sim21,sim22,sim23,sim24,sim25)
table2<-cbind(table2a[,1],table2a[,3],table2a[,2],table2a[,4],table2a[,5])
colnames(table2)<-c("BiasBint","CoverageBint","BiasBright","CoverageBright","Power")

write.csv(table2,"Simulation N=20000 12Pct Midpoint 7525 Results.csv")
write.csv(table2a,"Simulation N=20000 12Pct Midpoint 7525 Results Full.csv")


####################################################################################

set.seed(654)
nsims<-1000
sim1<-rowMeans(replicate(nsims,simRun2(n=20000,beta_right=-0.3,beta_int=-0.3,
	lambda_right=0.03,gammas_right=0.3,shape_int=2,scale_int=6)))
sim2<-rowMeans(replicate(nsims,simRun2(n=20000,beta_right=-0.3,beta_int=-0.15,
	lambda_right=0.03,gammas_right=0.3,shape_int=2,scale_int=6)))
sim3<-rowMeans(replicate(nsims,simRun2(n=20000,beta_right=-0.3,beta_int=0,
	lambda_right=0.03,gammas_right=0.3,shape_int=2,scale_int=6)))
sim4<-rowMeans(replicate(nsims,simRun2(n=20000,beta_right=-0.3,beta_int=0.15,
	lambda_right=0.03,gammas_right=0.3,shape_int=2,scale_int=6)))
sim5<-rowMeans(replicate(nsims,simRun2(n=20000,beta_right=-0.3,beta_int=0.3,
	lambda_right=0.03,gammas_right=0.3,shape_int=2,scale_int=6)))

sim6<-rowMeans(replicate(nsims,simRun2(n=20000,beta_right=-0.15,beta_int=-0.3,
	lambda_right=0.03,gammas_right=0.3,shape_int=2,scale_int=6)))
sim7<-rowMeans(replicate(nsims,simRun2(n=20000,beta_right=-0.15,beta_int=-0.15,
	lambda_right=0.03,gammas_right=0.3,shape_int=2,scale_int=6)))
sim8<-rowMeans(replicate(nsims,simRun2(n=20000,beta_right=-0.15,beta_int=0,
	lambda_right=0.03,gammas_right=0.3,shape_int=2,scale_int=6)))
sim9<-rowMeans(replicate(nsims,simRun2(n=20000,beta_right=-0.15,beta_int=0.15,
	lambda_right=0.03,gammas_right=0.3,shape_int=2,scale_int=6)))
sim10<-rowMeans(replicate(nsims,simRun2(n=20000,beta_right=-0.15,beta_int=0.3,
	lambda_right=0.03,gammas_right=0.3,shape_int=2,scale_int=6)))

sim11<-rowMeans(replicate(nsims,simRun2(n=20000,beta_right=0,beta_int=-0.3,
	lambda_right=0.03,gammas_right=0.3,shape_int=2,scale_int=6)))
sim12<-rowMeans(replicate(nsims,simRun2(n=20000,beta_right=0,beta_int=-0.15,
	lambda_right=0.03,gammas_right=0.3,shape_int=2,scale_int=6)))
sim13<-rowMeans(replicate(nsims,simRun2(n=20000,beta_right=0,beta_int=0,
	lambda_right=0.03,gammas_right=0.3,shape_int=2,scale_int=6)))
sim14<-rowMeans(replicate(nsims,simRun2(n=20000,beta_right=0,beta_int=0.15,
	lambda_right=0.03,gammas_right=0.3,shape_int=2,scale_int=6)))
sim15<-rowMeans(replicate(nsims,simRun2(n=20000,beta_right=0,beta_int=0.3,
	lambda_right=0.03,gammas_right=0.3,shape_int=2,scale_int=6)))

sim16<-rowMeans(replicate(nsims,simRun2(n=20000,beta_right=0.15,beta_int=-0.3,
	lambda_right=0.03,gammas_right=0.3,shape_int=2,scale_int=6)))
sim17<-rowMeans(replicate(nsims,simRun2(n=20000,beta_right=0.15,beta_int=-0.15,
	lambda_right=0.03,gammas_right=0.3,shape_int=2,scale_int=6)))
sim18<-rowMeans(replicate(nsims,simRun2(n=20000,beta_right=0.15,beta_int=0,
	lambda_right=0.03,gammas_right=0.3,shape_int=2,scale_int=6)))
sim19<-rowMeans(replicate(nsims,simRun2(n=20000,beta_right=0.15,beta_int=0.15,
	lambda_right=0.03,gammas_right=0.3,shape_int=2,scale_int=6)))
sim20<-rowMeans(replicate(nsims,simRun2(n=20000,beta_right=0.15,beta_int=0.3,
	lambda_right=0.03,gammas_right=0.3,shape_int=2,scale_int=6)))

sim21<-rowMeans(replicate(nsims,simRun2(n=20000,beta_right=0.3,beta_int=-0.3,
	lambda_right=0.03,gammas_right=0.3,shape_int=2,scale_int=6)))
sim22<-rowMeans(replicate(nsims,simRun2(n=20000,beta_right=0.3,beta_int=-0.15,
	lambda_right=0.03,gammas_right=0.3,shape_int=2,scale_int=6)))
sim23<-rowMeans(replicate(nsims,simRun2(n=20000,beta_right=0.3,beta_int=0,
	lambda_right=0.03,gammas_right=0.3,shape_int=2,scale_int=6)))
sim24<-rowMeans(replicate(nsims,simRun2(n=20000,beta_right=0.3,beta_int=0.15,
	lambda_right=0.03,gammas_right=0.3,shape_int=2,scale_int=6)))
sim25<-rowMeans(replicate(nsims,simRun2(n=20000,beta_right=0.3,beta_int=0.3,
	lambda_right=0.03,gammas_right=0.3,shape_int=2,scale_int=6)))

table2a<-rbind(sim1,sim2,sim3,sim4,sim5,sim6,sim7,sim8,sim9,sim10,sim11,sim12,sim13,sim14,sim15,
	sim16,sim17,sim18,sim19,sim20,sim21,sim22,sim23,sim24,sim25)
table2<-cbind(table2a[,1],table2a[,3],table2a[,2],table2a[,4],table2a[,5])
colnames(table2)<-c("BiasBint","CoverageBint","BiasBright","CoverageBright","Power")

write.csv(table2,"Simulation N=20000 12Pct Midpoint 2575 Results.csv")
write.csv(table2a,"Simulation N=20000 12Pct Midpoint 2575 Results Full.csv")


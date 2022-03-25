#make some plots to compare tables 1 and 2 for preventable interval censoring paper
#8/10/21

#comparing results for if you use midpt or upper value of interval within the composite outcome


#load libraries
library(blandr)
library(ggplot2)
library(gridExtra)




#7525 results from right and interval
upper.7525<-read.csv("Simulation N=20000 12Pct Upper 7525 Results.csv")
mdpt.7525<-read.csv("Simulation N=20000 12Pct Midpoint 7525 Results.csv")
plotbiasbint.7525<-blandr.draw(upper.7525$BiasBint,mdpt.7525$BiasBint)
plotcoveragebint.7525<-blandr.draw(upper.7525$CoverageBint,mdpt.7525$CoverageBint)
plotbiasbright.7525<-blandr.draw(upper.7525$BiasBright,mdpt.7525$BiasBright)
plotcoveragebright.7525<-blandr.draw(upper.7525$CoverageBright,mdpt.7525$CoverageBright)
plotpower.7525<-blandr.draw(upper.7525$Power,mdpt.7525$Power)



#5050 results 
upper.5050<-read.csv("Simulation N=20000 12Pct Upper 5050 Results.csv")
mdpt.5050<-read.csv("Simulation N=20000 12Pct Midpoint 5050 Results.csv")
plotbiasbint.5050<-blandr.draw(upper.5050$BiasBint,mdpt.5050$BiasBint)
plotcoveragebint.5050<-blandr.draw(upper.5050$CoverageBint,mdpt.5050$CoverageBint)
plotbiasbright.5050<-blandr.draw(upper.5050$BiasBright,mdpt.5050$BiasBright)
plotcoveragebright.5050<-blandr.draw(upper.5050$CoverageBright,mdpt.5050$CoverageBright)
plotpower.5050<-blandr.draw(upper.5050$Power,mdpt.5050$Power)



#2575 results from right and interval
upper.2575<-read.csv("Simulation N=20000 12Pct Upper 2575 Results.csv")
mdpt.2575<-read.csv("Simulation N=20000 12Pct Midpoint 2575 Results.csv")
plotbiasbint.2575<-blandr.draw(upper.2575$BiasBint,mdpt.2575$BiasBint)
plotcoveragebint.2575<-blandr.draw(upper.2575$CoverageBint,mdpt.2575$CoverageBint)
plotbiasbright.2575<-blandr.draw(upper.2575$BiasBright,mdpt.2575$BiasBright)
plotcoveragebright.2575<-blandr.draw(upper.2575$CoverageBright,mdpt.2575$CoverageBright)
plotpower.2575<-blandr.draw(upper.2575$Power,mdpt.2575$Power)


tiff("Bland Altman Plots Comparing Mdpt and Upper power.tiff",res=300,width=7,height=12,units="in")
grid.arrange(
	plotpower.5050+ggtitle("Power, 5050"),plotpower.7525+ggtitle("Power, 7525"),plotpower.2575+ggtitle("Power, 2575"),
	plotbiasbint.5050+ggtitle("Bias B Interval, 5050"),plotbiasbint.7525+ggtitle("Bias B Interval, 7525"),plotbiasbint.2575+ggtitle("Bias B Interval, 2575"),
	plotcoveragebint.5050+ggtitle("Coverage B Interval, 5050"),plotcoveragebint.7525+ggtitle("Coverage B Interval, 7525"),plotcoveragebint.2575+ggtitle("Coverage B Interval, 2575"),
	plotbiasbright.5050+ggtitle("Bias B Right, 5050"),plotbiasbright.7525+ggtitle("Bias B Right, 7525"),plotbiasbright.2575+ggtitle("Bias B Right, 2575"),
	plotcoveragebright.5050+ggtitle("Coverage B Right, 5050"),plotcoveragebright.7525+ggtitle("Coverage B Right, 7525"),plotcoveragebright.2575+ggtitle("Coverage B Right, 2575"),
	nrow=5,ncol=3)
dev.off()







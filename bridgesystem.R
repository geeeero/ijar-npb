# bridge system, 3 component types

#install.packages("devtools")
#library("devtools")
#install_github("louisaslett/ReliabilityTheory")
library("ReliabilityTheory")
library(ggplot2)
library(reshape2)

# produces survival signature matrix for one component of type "name",
# for use in nonParBayesSystemInference()
oneCompSurvSign <- function(name){
  res <- data.frame(name=c(0,1), Probability=c(0,1))
  names(res)[1] <- name
  res
}

# produces data frame with prior and posterior lower & upper component survival function
# for component of type "name" based on nonParBayesSystemInferencePriorSets() inputs
# for all components except survival signature; nLower, nUpper, yLower, yUpper must be data frames
# where each column corresponds to the component type, so there must be a match 
oneCompPriorPostSet <- function(name, at.times, test.data, nLower, nUpper, yLower, yUpper){
  sig <- oneCompSurvSign(name)
  nodata <- list(name=NULL)
  names(nodata) <- name
  nL <- nLower[, match(name, names(nLower))]
  nU <- nUpper[, match(name, names(nUpper))]
  yL <- yLower[, match(name, names(yLower))]
  yU <- yUpper[, match(name, names(yUpper))]
  data <- test.data[match(name, names(test.data))]
  prio <- nonParBayesSystemInferencePriorSets(at.times, sig, nodata, nL, nU, yL, yU)
  post <- nonParBayesSystemInferencePriorSets(at.times, sig,   data, nL, nU, yL, yU)
  data.frame(Time=rep(at.times,2), Lower=c(prio$lower,post$lower), Upper=c(prio$upper,post$upper),
             Item=rep(c("Prior", "Posterior"), each=length(at.times)))
}

# ----------------------------------------------

b3 <- graph.formula(s -- 1 -- 2:3 -- 4 -- 5:6 --t, 2 -- 5, 3 -- 6)
V(b3)$compType <- NA
V(b3)$compType[match(c("1"), V(b3)$name)] <- "T1"
V(b3)$compType[match(c("4"), V(b3)$name)] <- "T2"
V(b3)$compType[match(c("2", "3", "5", "6"), V(b3)$name)] <- "T3"
b3nulldata <- list("T1"=NULL, "T2"=NULL, "T3"=NULL)
b3testdata <- list("T1"=c(2.0, 2.2, 2.4, 2.6), "T2"=c(3.0, 3.2, 3.4, 3.6), "T3"=(1:4)/10+4) # T3 late failures
b3testdata <- list("T1"=c(2.0, 2.2, 2.4, 2.6), "T2"=c(3.0, 3.2, 3.4, 3.6), "T3"=(1:4)/10+0.5) # T3 early failures
b3testdata <- list("T1"=c(2.0, 2.2, 2.4, 2.6), "T2"=c(3.0, 3.2, 3.4, 3.6), "T3"=(1:4)-0.5) # T3 fitting failures
b3dat <- melt(b3testdata); names(b3dat) <- c("x", "Part")
b3sig <- computeSystemSurvivalSignature(b3)
b3t <- seq(0, 5, length.out=301)
b3nL <- data.frame(T1=rep(1,301), T2=rep(1,301), T3=rep(1,301))
b3nU <- data.frame(T1=rep(2,301), T2=rep(2,301), T3=rep(4,301))
b3yL <- data.frame(T1=rep(0.001, 301), T2=rep(0.001, 301), T3=c(rep(c(0.625,0.375,0.250,0.125,0.010), each=60), 0.01))
b3yU <- data.frame(T1=rep(0.999, 301), T2=rep(0.999, 301), T3=c(rep(c(0.990,0.875,0.500,0.375,0.250), each=60), 0.25))

b3T1 <- oneCompPriorPostSet("T1", b3t, b3testdata, b3nL, b3nU, b3yL, b3yU)
b3T2 <- oneCompPriorPostSet("T2", b3t, b3testdata, b3nL, b3nU, b3yL, b3yU)
b3T3 <- oneCompPriorPostSet("T3", b3t, b3testdata, b3nL, b3nU, b3yL, b3yU)
b3prio <- nonParBayesSystemInferencePriorSets(b3t, b3sig, b3nulldata, b3nL, b3nU, b3yL, b3yU)
b3post <- nonParBayesSystemInferencePriorSets(b3t, b3sig, b3testdata, b3nL, b3nU, b3yL, b3yU)

b3df <- rbind(data.frame(b3T1, Part="T1"), data.frame(b3T2, Part="T2"), data.frame(b3T3, Part="T3"),
              data.frame(Time=rep(b3t,2), Lower=c(b3prio$lower,b3post$lower), Upper=c(b3prio$upper,b3post$upper),
                         Item=rep(c("Prior", "Posterior"), each=length(b3t)), Part="System"))
b3df$Item <- ordered(b3df$Item, levels=c("Prior", "Posterior"))

p3 <- ggplot(b3df) + geom_ribbon(aes(x=Time, ymin=Lower, ymax=Upper, group=Item, colour=Item, fill=Item), alpha=0.3)
p3 <- p3 + facet_wrap(~Part, ncol=2) + geom_rug(aes(x=x), data=g3dat) + xlab("Time") + ylab("Survival Probability")
#pdf("3comp-latefailures.pdf", width=8, height=5)
#pdf("3comp-earlyfailures.pdf", width=8, height=5)
pdf("3comp-fittingfailures.pdf", width=8, height=5)
p3
dev.off()

#
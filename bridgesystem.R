# bridge system, 3 component types

#install.packages("devtools")
#library("devtools")
#install_github("louisaslett/ReliabilityTheory")
library("ReliabilityTheory")
library(ggplot2)
library(reshape2)
library(RColorBrewer)
library(xtable)

bottomlegend <- theme(legend.position = 'bottom', legend.direction = 'horizontal', legend.title = element_blank())
rightlegend <- theme(legend.title = element_blank())

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

b3 <- graph.formula(s -- 2:3 -- 4 -- 5:6 -- 1 -- t, 2 -- 5, 3 -- 6)
b3 <- setCompTypes(b3, list("T1"=c(2,3,5,6), "T2"=c(4), "T3"=c(1)))
b3nulldata <- list("T1"=NULL, "T2"=NULL, "T3"=NULL)
b3testdata <- list("T1"=c(2.2, 2.4, 2.6, 2.8), "T2"=c(3.2, 3.4, 3.6, 3.8), "T3"=(1:4)/10+4) # T3 late failures
b3testdata <- list("T1"=c(2.2, 2.4, 2.6, 2.8), "T2"=c(3.2, 3.4, 3.6, 3.8), "T3"=(1:4)/10+0.5) # T3 early failures
b3testdata <- list("T1"=c(2.2, 2.4, 2.6, 2.8), "T2"=c(3.2, 3.4, 3.6, 3.8), "T3"=(1:4)-0.5) # T3 fitting failures
b3dat <- melt(b3testdata); names(b3dat) <- c("x", "Part")
b3dat$Part <- ordered(b3dat$Part, levels=c("T1", "T2", "T3", "System"))
b3sig <- computeSystemSurvivalSignature(b3)
b3t <- seq(0, 5, length.out=301)
b3nL <- data.frame(T1=rep(1,301), T2=rep(1,301), T3=rep(1,301))
b3nU <- data.frame(T1=rep(2,301), T2=rep(2,301), T3=rep(4,301))
b3yL <- data.frame(T1=rep(0.001, 301), T2=rep(0.001, 301), T3=c(rep(c(0.625,0.375,0.250,0.125,0.010), each=60), 0.01))
b3yU <- data.frame(T1=rep(0.999, 301), T2=rep(0.999, 301), T3=c(rep(c(0.999,0.875,0.500,0.375,0.250), each=60), 0.25))

b3T1 <- oneCompPriorPostSet("T1", b3t, b3testdata, b3nL, b3nU, b3yL, b3yU)
b3T2 <- oneCompPriorPostSet("T2", b3t, b3testdata, b3nL, b3nU, b3yL, b3yU)
b3T3 <- oneCompPriorPostSet("T3", b3t, b3testdata, b3nL, b3nU, b3yL, b3yU)
b3prio <- nonParBayesSystemInferencePriorSets(b3t, b3sig, b3nulldata, b3nL, b3nU, b3yL, b3yU)
b3post <- nonParBayesSystemInferencePriorSets(b3t, b3sig, b3testdata, b3nL, b3nU, b3yL, b3yU)

b3df <- rbind(data.frame(b3T1, Part="T1"), data.frame(b3T2, Part="T2"), data.frame(b3T3, Part="T3"),
              data.frame(Time=rep(b3t,2), Lower=c(b3prio$lower,b3post$lower), Upper=c(b3prio$upper,b3post$upper),
                         Item=rep(c("Prior", "Posterior"), each=length(b3t)), Part="System"))
b3df$Item <- ordered(b3df$Item, levels=c("Prior", "Posterior"))
b3df$Part <- ordered(b3df$Part, levels=c("T1", "T2", "T3", "System"))

# + scale_colour_manual(values = c("red","blue"))
# + scale_fill_brewer(palette="Set3")
p3 <- ggplot(b3df, aes(x=Time)) + theme_bw()
p3 <- p3 + scale_fill_manual(values = c("#b2df8a", "#1f78b4")) + scale_colour_manual(values = c("#b2df8a", "#1f78b4"))
p3 <- p3 + geom_line(aes(y=Upper, group=Item, colour=Item)) + geom_line(aes(y=Lower, group=Item, colour=Item))
p3 <- p3 + geom_ribbon(aes(ymin=Lower, ymax=Upper, group=Item, colour=Item, fill=Item), alpha=0.5)
p3 <- p3  + facet_wrap(~Part, nrow=2) + geom_rug(aes(x=x), data=b3dat) + xlab("Time") + ylab("Survival Probability")
pdf("bridge-latefailures.pdf", width=8, height=5)
pdf("bridge-earlyfailures.pdf", width=8, height=5)
pdf("bridge-fittingfailures.pdf", width=8, height=5)
p3 + rightlegend
dev.off()

b3sigtable <- b3sig[b3sig$T3 == 1,]
b3sigtable$T1 <- as.factor(b3sigtable$T1)
b3sigtable$T2 <- as.factor(b3sigtable$T2)
b3sigtable$T3 <- as.factor(b3sigtable$T3)
xtable(b3sigtable)

#
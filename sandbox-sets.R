#play around with sets

#install.packages("devtools")
#library("devtools")
#install_github("louisaslett/ReliabilityTheory")
library("ReliabilityTheory")
library(ggplot2)
library(reshape2)


# Risk Analysis paper example with sets of priors
# First specify the system layout, numbered as per Figure 4
g <- graph.formula(s -- 1 -- 2:4:5, 2 -- 3 -- t, 4:5 -- 6 -- t,
                   s -- 7 -- 8 -- t, s -- 9 -- 10 -- 11 -- t, 7 -- 10 -- 8)
# Then specify the component types (within the circles of Figure 4)
V(g)$compType <- NA
V(g)$compType[match(c("1","6","11"), V(g)$name)] <- "T1"
V(g)$compType[match(c("2","3","9"), V(g)$name)] <- "T2"
V(g)$compType[match(c("4","5","10"), V(g)$name)] <- "T3"
V(g)$compType[match(c("7","8"), V(g)$name)] <- "T4"
# Compute the survival signature table from Appendix
sig <- computeSystemSurvivalSignature(g)
# Simulate the test data (same seed as used in the paper)
set.seed(233)
t1 <- rexp(100, rate=0.55)
t2 <- rweibull(100, scale=1.8, shape=2.2)
t3 <- rlnorm(100, 0.4, 0.9)
t4 <- rgamma(100, scale=0.9, shape=3.2)
# Compile into a list as required by this function
test.data <- list("T1"=t1, "T2"=t2, "T3"=t3, "T4"=t4)
no.test.data <- list("T1"=NULL, "T2"=NULL, "T3"=NULL, "T4"=NULL)
# Create a vector of times at which to evaluate the posterior predictive
# survival probability and compute using this function
t <- seq(0, 5, length.out=300)
yS <- nonParBayesSystemInference(t, sig, test.data) # with fixed prior as before
# Set of prior R(t)'s
yS0 <- nonParBayesSystemInferencePriorSets(t, sig, no.test.data, nLower=1, nUpper=50, yLower=0.4, yUpper=0.7)
# Set of posterior R(t)'s
yS2 <- nonParBayesSystemInferencePriorSets(t, sig, test.data, nLower=1, nUpper=50, yLower=0.4, yUpper=0.7)
# Plot
p <- ggplot(data.frame(Time=rep(t,3), Probability=c(yS, yS2$lower, yS2$upper),
                       Item=c(rep(c("System", "Lower", "Upper"), each=300))))
p <- p + geom_line(aes(x=Time, y=Probability, linetype=Item))
p <- p + xlab("Time") + ylab("Survival Probability")
p

p2 <- ggplot(data.frame(Time=rep(t,2), Lower=c(yS0$lower,yS2$lower), Upper=c(yS0$upper,yS2$upper),
                        Item=rep(c("Prior", "Posterior"), each=300)))
p2 <- p2 + geom_ribbon(aes(x=Time, ymin=Lower, ymax=Upper, group=Item, colour=Item, fill=Item), alpha=0.3)
p2 <- p2 + xlab("Time") + ylab("Survival Probability")
p2

# ----------------------------------------------

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

# parallel system with two components, T1 prior informative, T2 prior near ignorance
g2p <- graph.formula(s -- 1:2 -- t)
V(g2p)$compType <- NA
V(g2p)$compType[match(c("1"), V(g2p)$name)] <- "T1"
V(g2p)$compType[match(c("2"), V(g2p)$name)] <- "T2"
g2pnulldata <- list("T1"=NULL, "T2"=NULL)
g2ptestdata <- list("T1"=1:4, "T2"=1:4+0.5)
g2pdat <- rbind(data.frame(x=g2ptestdata$"T1", Part="T1"), data.frame(x=g2ptestdata$"T2", Part="T2"))
g2psig <- computeSystemSurvivalSignature(g2p)
g2pt <- seq(0, 5, length.out=11)
g2nL <- data.frame(T1=rep(1,11), T2=rep(1,11))
g2nU <- data.frame(T1=rep(10,11), T2=rep(2,11))
g2yL <- data.frame(T1=c(0.75, rep(c(0.75,0.50,0.01), each=3), 0.01), T2=rep(0.01, 11))
g2yU <- data.frame(T1=c(0.99, rep(c(0.99,0.75,0.50), each=3), 0.50), T2=rep(0.99, 11))

g2sT1 <- oneCompPriorPostSet("T1", g2pt, g2ptestdata, g2nL, g2nU, g2yL, g2yU)
g2sT2 <- oneCompPriorPostSet("T2", g2pt, g2ptestdata, g2nL, g2nU, g2yL, g2yU)
#ggplot(g2sT1) + geom_ribbon(aes(x=Time, ymin=Lower, ymax=Upper, group=Item, colour=Item, fill=Item), alpha=0.3)
#ggplot(g2sT2) + geom_ribbon(aes(x=Time, ymin=Lower, ymax=Upper, group=Item, colour=Item, fill=Item), alpha=0.3)
g2sprio <- nonParBayesSystemInferencePriorSets(g2pt, g2psig, g2pnulldata, g2nL, g2nU, g2yL, g2yU)
g2spost <- nonParBayesSystemInferencePriorSets(g2pt, g2psig, g2ptestdata, g2nL, g2nU, g2yL, g2yU)

g2sdf <- rbind(data.frame(g2sT1, Part="T1"), data.frame(g2sT2, Part="T2"),
               data.frame(Time=rep(g2pt,2), Lower=c(g2sprio$lower,g2spost$lower), Upper=c(g2sprio$upper,g2spost$upper),
                          Item=rep(c("Prior", "Posterior"), each=length(g2pt)), Part="System"))

p2s <- ggplot(g2sdf) + geom_ribbon(aes(x=Time, ymin=Lower, ymax=Upper, group=Item, colour=Item, fill=Item), alpha=0.3)
p2s <- p2s + facet_wrap(~Part, ncol=2) + geom_rug(aes(x=x), data=g2pdat) + xlab("Time") + ylab("Survival Probability")
p2s

# ----------------------------------------------

# 2-parallel system with each subsystem 3 component series system
g3 <- graph.formula(s -- 1 -- 2 -- 3 -- t, s -- 4 -- 5 -- 6 -- t)
V(g3)$compType <- NA
V(g3)$compType[match(c("1", "4"), V(g3)$name)] <- "T1"
V(g3)$compType[match(c("2", "5"), V(g3)$name)] <- "T2"
V(g3)$compType[match(c("3", "6"), V(g3)$name)] <- "T3"
g3nulldata <- list("T1"=NULL, "T2"=NULL, "T3"=NULL)
g3testdata <- list("T1"=c(2, 2.2, 2.4, 2.6), "T2"=c(3.0, 3.2, 3.4, 3.6), "T3"=(1:4)/10+4) # T3 late failures
g3testdata <- list("T1"=c(2, 2.2, 2.4, 2.6), "T2"=c(3.0, 3.2, 3.4, 3.6), "T3"=(1:4)/10+0.5) # T3 early failures
g3testdata <- list("T1"=c(2, 2.2, 2.4, 2.6), "T2"=c(3.0, 3.2, 3.4, 3.6), "T3"=(1:4)-0.5) # T3 fitting failures
g3dat <- melt(g3testdata); names(g3dat) <- c("x", "Part")
g3sig <- computeSystemSurvivalSignature(g3)
g3t <- seq(0, 5, length.out=301)
g3nL <- data.frame(T1=rep(1,301), T2=rep(1,301), T3=rep(1,301))
g3nU <- data.frame(T1=rep(2,301), T2=rep(2,301), T3=rep(4,301))
g3yL <- data.frame(T1=rep(0.01, 301), T2=rep(0.01, 301), T3=c(rep(c(0.625,0.375,0.250,0.125,0.010), each=60), 0.01))
g3yU <- data.frame(T1=rep(0.99, 301), T2=rep(0.99, 301), T3=c(rep(c(0.990,0.875,0.500,0.375,0.250), each=60), 0.25))

g3T1 <- oneCompPriorPostSet("T1", g3t, g3testdata, g3nL, g3nU, g3yL, g3yU)
g3T2 <- oneCompPriorPostSet("T2", g3t, g3testdata, g3nL, g3nU, g3yL, g3yU)
g3T3 <- oneCompPriorPostSet("T3", g3t, g3testdata, g3nL, g3nU, g3yL, g3yU)
g3prio <- nonParBayesSystemInferencePriorSets(g3t, g3sig, g3nulldata, g3nL, g3nU, g3yL, g3yU)
g3post <- nonParBayesSystemInferencePriorSets(g3t, g3sig, g3testdata, g3nL, g3nU, g3yL, g3yU)

g3df <- rbind(data.frame(g3T1, Part="T1"), data.frame(g3T2, Part="T2"), data.frame(g3T3, Part="T3"),
              data.frame(Time=rep(g3t,2), Lower=c(g3prio$lower,g3post$lower), Upper=c(g3prio$upper,g3post$upper),
                         Item=rep(c("Prior", "Posterior"), each=length(g3t)), Part="System"))
g3df$Item <- ordered(g3df$Item, levels=c("Prior", "Posterior"))

p3 <- ggplot(g3df) + geom_ribbon(aes(x=Time, ymin=Lower, ymax=Upper, group=Item, colour=Item, fill=Item), alpha=0.3)
p3 <- p3 + facet_wrap(~Part, ncol=2) + geom_rug(aes(x=x), data=g3dat) + xlab("Time") + ylab("Survival Probability")
#pdf("3comp-latefailures.pdf", width=8, height=5)
#pdf("3comp-earlyfailures.pdf", width=8, height=5)
pdf("3comp-fittingfailures.pdf", width=8, height=5)
p3
dev.off()

# ----------------------------------------------

#
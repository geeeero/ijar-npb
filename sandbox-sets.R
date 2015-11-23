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
oneCompPriorPost <- function(name, at.times, test.data, alpha, beta){
  sig <- oneCompSurvSign(name)
  nodata <- list(name=NULL)
  names(nodata) <- name
  avec <- alpha[, match(name, names(alpha))]
  bvec <- beta[, match(name, names(beta))]
  data <- test.data[match(name, names(test.data))]
  prio <- nonParBayesSystemInference(at.times, sig, nodata, avec, bvec)
  post <- nonParBayesSystemInference(at.times, sig, data, avec, bvec)
  data.frame(Time=at.times, Prior=prio, Posterior=post)
}

#
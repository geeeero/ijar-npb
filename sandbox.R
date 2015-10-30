# playing around

setwd("../ijar-npb")
library("ReliabilityTheory")
library(ggplot2)

# Risk Analysis paper example
g <- graph.formula(s -- 1 -- 2:4:5, 2 -- 3 -- t, 4:5 -- 6 -- t,
                   s -- 7 -- 8 -- t, s -- 9 -- 10 -- 11 -- t, 7 -- 10 -- 8)
V(g)$compType <- NA
V(g)$compType[match(c("1","6","11"), V(g)$name)] <- "T1"
V(g)$compType[match(c("2","3","9"), V(g)$name)] <- "T2"
V(g)$compType[match(c("4","5","10"), V(g)$name)] <- "T3"
V(g)$compType[match(c("7","8"), V(g)$name)] <- "T4"

sig <- computeSystemSurvivalSignature(g)

set.seed(233)
t1 <- rexp(100, rate=0.55)
t2 <- rweibull(100, scale=1.8, shape=2.2)
t3 <- rlnorm(100, 0.4, 0.9)
t4 <- rgamma(100, scale=0.9, shape=3.2)

d <- ggplot(data.frame(Failuretimes=c(t1,t2,t3,t4),
                       Item=c(rep(c("T1", "T2", "T3", "T4"), each=100))))
d <- d + geom_histogram(aes(x=Failuretimes), binwidth=0.5) + facet_grid(Item ~ .)
d <- d + coord_cartesian(xlim=c(0,11)) + xlab("Failure times")
d

t <- seq(0, 5, length.out=300)

# prior survival function
no.test.data <- list("T1"=NULL, "T2"=NULL, "T3"=NULL, "T4"=NULL)
s0 <- nonParBayesSystemInference(t, sig, no.test.data)
pr <- ggplot(data.frame(Time=t, Probability=s0))
pr <- pr + geom_line(aes(x=Time, y=Probability)) + coord_cartesian(ylim=c(0,1))
pr <- pr + xlab("Time") + ylab("Survival Probability")
pr

# posterior survival function
test.data <- list("T1"=t1, "T2"=t2, "T3"=t3, "T4"=t4)
yS <- nonParBayesSystemInference(t, sig, test.data)
p <- ggplot(data.frame(Time=t, Probability=yS))
p <- p + geom_line(aes(x=Time, y=Probability))
p <- p + xlab("Time") + ylab("Survival Probability")
p


# ----------------------------------------------

g <- graph.formula(s -- 1:2:3 -- t)
V(g)$compType <- NA
V(g)$compType[match(c("1"), V(g)$name)] <- "T1"
V(g)$compType[match(c("2", "3"), V(g)$name)] <- "T2"
no.test.data <- list("T1"=NULL, "T2"=NULL)
test.data <- list("T1"=t1, "T2"=t2)
sig <- computeSystemSurvivalSignature(g)
t <- seq(0, 5, length.out=300)
s0 <- nonParBayesSystemInference(t, sig, no.test.data)
yS <- nonParBayesSystemInference(t, sig, test.data)

pplot <- function() {
  p <- ggplot(data.frame(Time=t, Prior=s0, Posterior=yS))
  p <- p + geom_line(aes(x=Time, y=Prior)) + geom_line(aes(x=Time, y=Posterior))
  p <- p + xlab("Time") + ylab("Survival Probability") + coord_cartesian(ylim=c(0,1))
  p
}
pplot()

g <- graph.formula(s -- 1 -- t)
V(g)$compType <- NA
V(g)$compType[match(c("1"), V(g)$name)] <- "T1"
sig <- computeSystemSurvivalSignature(g)
no.test.data <- list("T1"=NULL)
test.data <- list("T1"=1:4)
t <- seq(0, 5, length.out=11)
s0 <- nonParBayesSystemInference(t, sig, no.test.data) # error fixed
yS <- nonParBayesSystemInference(t, sig, test.data)    # error fixed
pplot()

avec <- c(rep(5:1, each=2), 1)
bvec <- c(rep(5:1, each=2), 1)
s0 <- nonParBayesSystemInference(t, sig, no.test.data, avec, bvec)
yS <- nonParBayesSystemInference(t, sig, test.data, avec, bvec)
pplot()


#

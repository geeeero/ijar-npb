# playing around

#setwd("../ijar-npb")
#install.packages("devtools")
#library("devtools")
#install_github("louisaslett/ReliabilityTheory")
library("ReliabilityTheory")
library(ggplot2)
library(reshape2)

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
  p <- p + xlab("Time") + ylab("Survival Probability") + coord_cartesian(ylim=c(-0.05,1.05))
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

# ----------------------------------------------

# converts n, y parameters to Beta parameters alpha, beta
ny2ab <- function(n,y){
  a <- n*y
  b <- n*(1-y)
  data.frame(a=a, b=b)
}

# produces survival signature matrix for one component of type "name",
# for use in nonParBayesSystemInference()
oneCompSurvSign <- function(name){
  res <- data.frame(name=c(0,1), Probability=c(0,1))
  names(res)[1] <- name
  res
}

# produces data frame with prior and posterior component survival function
# for component of type "name" based on nonParBayesSystemInference() inputs
# for all components (except survival signature; alpha and beta must be data frames)
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



# ----------------------------------------------

g2p <- graph.formula(s -- 1:2 -- t)
V(g2p)$compType <- NA
V(g2p)$compType[match(c("1"), V(g2p)$name)] <- "T1"
V(g2p)$compType[match(c("2"), V(g2p)$name)] <- "T2"
g2pnulldata <- list("T1"=NULL, "T2"=NULL)
g2ptestdata <- list("T1"=1:4, "T2"=1:4)
g2psig <- computeSystemSurvivalSignature(g2p)
g2pt <- seq(0, 5, length.out=11)
g2T1ab <- ny2ab(rep(1,11), seq(0.99, 0.01, length.out=11))
g2adf <- data.frame("T1"=g2T1ab$a, "T2"=rep(1,11))
g2bdf <- data.frame("T1"=g2T1ab$b, "T2"=rep(1,11))

g2pprio <- nonParBayesSystemInference(g2pt, g2psig, g2pnulldata, g2adf, g2bdf)
g2ppost <- nonParBayesSystemInference(g2pt, g2psig, g2ptestdata, g2adf, g2bdf)

g2 <- ggplot(data.frame(Time=rep(g2pt, 2), Surv=c(g2pprio, g2ppost),
                        What=rep(c("Prior", "Posterior"), each=length(g2pt))))
g2 <- g2 + geom_line(aes(x=Time, y=Surv, linetype=What))
g2 <- g2 + ylab("Survival Probability") + coord_cartesian(ylim=c(-0.05,1.05))
g2

# look at component priors and posteriors
g2T1 <- oneCompPriorPost("T1", g2pt, g2ptestdata, g2adf, g2bdf)
ggplot(g2T1) + geom_line(aes(x=Time, y=Prior)) + geom_line(aes(x=Time, y=Posterior))
g2T2 <- oneCompPriorPost("T2", g2pt, g2ptestdata, g2adf, g2bdf)
ggplot(g2T2) + geom_line(aes(x=Time, y=Prior)) + geom_line(aes(x=Time, y=Posterior))

g2df <- data.frame(Time=g2pt, SysPrior=g2pprio, SysPosterior=g2ppost,
                   T1Prior=g2T1$Prior, T1Posterior=g2T1$Posterior,
                   T2Prior=g2T2$Prior, T2Posterior=g2T2$Posterior)
g2df <- melt(g2df, "Time")
g2 <- ggplot(g2df) + geom_line(aes(x=Time, y=value, color=variable))
g2 <- g2 + ylab("Survival Probability") + coord_cartesian(ylim=c(-0.05,1.05))
g2

bottomlegend <- theme(legend.position = 'bottom', legend.direction = 'horizontal', legend.title = element_blank())
rightlegend <- theme(legend.title = element_blank())

g2df <- rbind(cbind(melt(g2T1, "Time"), Part="T1"), cbind(melt(g2T2, "Time"), Part="T2"),
              cbind(melt(data.frame(Time=g2pt, Prior=g2pprio, Posterior=g2ppost), "Time"), Part="System"))
g2 <- ggplot(g2df) + geom_line(aes(x=Time, y=value, color=Part, linetype=variable), size=1)
g2 <- g2 + ylab("Survival Probability") + # coord_cartesian(ylim=c(-0.05,1.05)) +
  rightlegend + # theme_minimal() +
#  scale_colour_manual(values=c("red", "orange", "yellow"))
#  scale_colour_brewer(palette="Set1")
  #  scale_colour_manual(values=rainbow(3))
  scale_colour_manual(values=heat.colors(3))
g2

g2 <- ggplot(g2df) + geom_line(aes(x=Time, y=value, linetype=variable)) + 
  geom_rug(aes(x=x), data=rbind(data.frame(x=g2ptestdata$"T1", Part="T1"), data.frame(x=g2ptestdata$"T2", Part="T2"))) +
#  facet_grid(Part ~ .) # arranges facets in a column
  facet_wrap(~ Part, ncol=2) # for, e.g., 3 component types, this will arrange the 4 facets in 2x2 
g2 <- g2 + ylab("Survival Probability") + bottomlegend
g2

# ----------------------------------------------


#

# small automobile breaking system example form Ebeling (1997, pp. 91f)

#install.packages("devtools")
#library("devtools")
#install_github("louisaslett/ReliabilityTheory")
library("ReliabilityTheory")
library(ggplot2)
library(reshape2)

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


ab <- graph.formula(s -- M -- C1:C2:C3:C4, P1:P2:P3:P4 -- t,
                    C1 -- P1, C2 -- P2, C3 -- P3, C4 -- P4, s -- H -- P3:P4)
ab <- setCompTypes(ab, list("M"=c("M"), "H"=c("H"), "C"=c("C1", "C2", "C3", "C4"), "P"=c("P1", "P2", "P3", "P4")))
# data
set.seed(233)
Mdata <- rexp(5, rate=0.25)
Hdata <- rlnorm(10, 1.5, 0.3)
Cdata <- rexp(15, rate=0.3)
Pdata <- rgamma(20, scale=0.9, shape=3.2)
abnulldata <- list("M"=NULL, "H"=NULL, "C"=NULL, "P"=NULL)
abtestdata <- list("M"=Mdata, "H"=Hdata, "C"=Cdata, "P"=Pdata)
abdat <- melt(abtestdata); names(abdat) <- c("x", "Part")
absig <- computeSystemSurvivalSignature(ab)
abt <- seq(0, 10, length.out=301)
# priors
abnL <- data.frame(M=rep(1,301), H=rep(1,301), C=rep(1,301), P=rep(1,301))
abnU <- data.frame(M=rep(8,301), H=rep(2,301), C=rep(2,301), P=rep(2,301))
abyL <- data.frame(M=c(rep(0.8, 150), rep(0.6, 60), rep(0.2, 30), rep(0.1, 61)),
                   H=c(rep(0.5, 150), rep(0.25, 60), rep(0.01, 91)),
                   C=c(rep(c(0.75, 0.73, 0.71, 0.70, 0.60, 0.45, 0.30, 0.23, 0.21, 0.20), each=30), 0.20),
                   P=c(rep(c(0.5, 0.01), each=150), 0.01))
abyU <- data.frame(M=c(rep(0.99, 180), rep(0.9, 60), rep(0.6, 30), rep(0.4, 31)),
                   H=c(rep(0.99, 90), rep(0.9, 90), rep(0.7, 30), rep(0.5, 30), rep(0.3,61)),
                   C=c(rep(c(0.99, 0.98, 0.96, 0.95, 0.90, 0.65, 0.50, 0.45, 0.43, 0.42), each=30), 0.42),
                   P=c(rep(c(0.99, 0.65), each=150), 0.65))
#posteriors
abM <- oneCompPriorPostSet("M", abt, abtestdata, abnL, abnU, abyL, abyU)
abH <- oneCompPriorPostSet("H", abt, abtestdata, abnL, abnU, abyL, abyU)
abC <- oneCompPriorPostSet("C", abt, abtestdata, abnL, abnU, abyL, abyU)
abP <- oneCompPriorPostSet("P", abt, abtestdata, abnL, abnU, abyL, abyU)
abprio <- nonParBayesSystemInferencePriorSets(abt, absig, abnulldata, abnL, abnU, abyL, abyU)
abpost <- nonParBayesSystemInferencePriorSets(abt, absig, abtestdata, abnL, abnU, abyL, abyU)
#data frame for plot
abdf <- rbind(data.frame(abM, Part="M"), data.frame(abH, Part="H"), data.frame(abC, Part="C"), data.frame(abP, Part="P"),
              data.frame(Time=rep(abt,2), Lower=c(abprio$lower,abpost$lower), Upper=c(abprio$upper,abpost$upper),
                         Item=rep(c("Prior", "Posterior"), each=length(abt)), Part="System"))
abdf$Item <- ordered(abdf$Item, levels=c("Prior", "Posterior"))
#the plot
ab1 <- ggplot(abdf, aes(x=Time))
ab1 <- ab1 +  geom_line(aes(y=Upper, group=Item, colour=Item)) + geom_line(aes(y=Lower, group=Item, colour=Item))
ab1 <- ab1 + geom_ribbon(aes(ymin=Lower, ymax=Upper, group=Item, colour=Item, fill=Item), alpha=0.3)
ab1 <- ab1 + facet_wrap(~Part, nrow=2) + geom_rug(aes(x=x), data=abdat) + xlab("Time") + ylab("Survival Probability")
ab1 <- ab1 + bottomlegend

pdf("breakingsystem-2.pdf", width=8, height=5)
ab1
dev.off()
#
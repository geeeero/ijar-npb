# Illustration prior-data conflict in Beta-Binomial model
# plots of posterior predictive

library(ggplot2)
library(reshape2)
library(grid)
library(gridExtra)
library(luck)
source("../../lund_1512/lund-1512/course/04-01_BinomialData.R")
source("../../lund_1512/lund-1512/course/04-02_Binomial.R")

updateLuckY <- function (n0, y0, tau, n){ (n0*y0+tau)/(n0+n) }
updateLuckN <- function (n0, n){ n0+n }

nyupdate <- function (pr, data){
  nn <- updateLuckN(pr[1], data[2])
  yn <- updateLuckY(pr[1], pr[2], data[1], data[2])
  c(nn,yn)
}

luck4cny <- function(luck, posterior=FALSE){
  c1 <- c(n0(luck)[1], y0(luck)[1])
  c2 <- c(n0(luck)[1], y0(luck)[2])
  c3 <- c(n0(luck)[2], y0(luck)[1])
  c4 <- c(n0(luck)[2], y0(luck)[2])
  if(posterior){
    c1 <- nyupdate(c1, c(tau(data(luck)), n(data(luck))))
    c2 <- nyupdate(c2, c(tau(data(luck)), n(data(luck))))
    c3 <- nyupdate(c3, c(tau(data(luck)), n(data(luck))))
    c4 <- nyupdate(c4, c(tau(data(luck)), n(data(luck))))
  }
  list(c1=c1, c2=c2, c3=c3, c4=c4)
}

dbetany <- function(x, ny, ...){
  dbeta(x, shape1=ny[1]*ny[2], shape2=ny[1]*(1-ny[2]), ...)
}

pbetany <- function(x, ny, ...){
  pbeta(x, shape1=ny[1]*ny[2], shape2=ny[1]*(1-ny[2]), ...)
}

dbebinvec <- function(pos, m, l = NULL){
  if (is.null(l))
    l <- 0:m
  sapply(l, function(li, m, n, y){
    choose(m,li) * beta(li+n*y, m-li+n*(1-y)) / beta(n*y, n*(1-y)) 
  }, m=m, n=pos[1], y=pos[2])
}

pbebinvec <- function(pos, m, l = NULL){
  res <- dbebinvec(pos=pos, m=m, l=NULL)
  res <- cumsum(res)
  if(!is.null(l))
    return(res[l+1])
  else
    return(res)
}

# one prior updated to posteriors with same variance
sgpr <- c(8, 0.75) # n0, y0
data1 <- c(12,16)  #  s,  n
data2 <- c( 0,16)

pos1 <- nyupdate(sgpr, data1)
pos2 <- nyupdate(sgpr, data2)

betavec <- seq(0,1, length.out=200)
# beta densities
bdsgpr <- dbetany(betavec, sgpr)
bdpos1 <- dbetany(betavec, pos1)
bdpos2 <- dbetany(betavec, pos2)
# beta cdfs
bpsgpr <- pbetany(betavec, sgpr)
bppos1 <- pbetany(betavec, pos1)
bppos2 <- pbetany(betavec, pos2)

bottomlegend <- theme(legend.position = 'bottom', legend.direction = 'horizontal', legend.title = element_blank())
rightlegend <- theme(legend.title = element_blank())

bddf <- data.frame(x=betavec, Prior=bdsgpr, "Posterior1"=bdpos1, "Posterior2"=bdpos2)
bpdf <- data.frame(x=betavec, Prior=bpsgpr, "Posterior1"=bppos1, "Posterior2"=bppos2)
names(bddf)[c(3,4)] <- c("Posterior 1", "Posterior 2")
names(bpdf)[c(3,4)] <- c("Posterior 1", "Posterior 2")
bddf <- melt(bddf, "x")
bpdf <- melt(bpdf, "x")
bdf <- rbind(data.frame(bddf, Item="pdf"), data.frame(bpdf, Item="cdf"))

betadens <- ggplot(bddf, aes(x, value, group=variable, linetype=variable)) + geom_line() +
  bottomlegend + xlab(expression(p[t]^k)) + ylab(expression(f(p[t]^k)))
betadens

betadens2 <- ggplot(bddf, aes(x, value, group=variable, linetype=variable)) + geom_line() +
  guides(linetype="none") + xlab(expression(p[t]^k)) + ylab(expression(f(p[t]^k)))
betadens2

betapdfcdf <- ggplot(bdf, aes(x, value, group=variable, linetype=variable)) + geom_line() +
  guides(linetype="none") + xlab(expression(p[t]^k)) + ylab("") +
  facet_grid(Item ~ ., scales="free")
betapdfcdf


# beta-binomial distributions
#dbebinvec(pos1,2,1)
m <- 5
bebinpr <- dbebinvec(sgpr, m)
bebin1 <- dbebinvec(pos1, m)
bebin2 <- dbebinvec(pos2, m)
bebinprc <- cumsum(bebinpr)
bebin1c <- cumsum(bebin1)
bebin2c <- cumsum(bebin2)

pdcilludf <- data.frame(x=0:m, Prior=bebinpr, "Posterior1"=bebin1, "Posterior2"=bebin2)
names(pdcilludf)[c(3,4)] <- c("Posterior 1", "Posterior 2")
pdcilludf$x <- ordered(pdcilludf$x)
pdcilludf <- melt(pdcilludf, "x")
pdcillu <- ggplot(pdcilludf, aes(x, value, group=variable, linetype=variable)) + geom_point() + geom_line() +
  rightlegend + xlab(expression(l[k])) + ylab(expression(P(C[t]^k == l[k])))
pdcillu

pdcilludfc <- data.frame(x=0:m, Prior=bebinprc, "Posterior1"=bebin1c, "Posterior2"=bebin2c)
names(pdcilludfc)[c(3,4)] <- c("Posterior 1", "Posterior 2")
pdcilludfc$x <- ordered(pdcilludfc$x)
pdcilludfc <- melt(pdcilludfc, "x")

pdcilludfdfc <- rbind(data.frame(pdcilludf, Item="pmf"), data.frame(pdcilludfc, Item="cmf"))

pdcillupmfcmf <- ggplot(pdcilludfdfc, aes(x, value, group=variable, linetype=variable)) + geom_point() + geom_line() +
  rightlegend + xlab(expression(l[k])) + ylab("") + facet_grid(Item ~ ., scales="free")
pdcillupmfcmf


pdcilludf2 <- data.frame(lk=rep(0:m, 3), plk=c(bebinpr, bebin1, bebin2),
                        which=rep(c("Prior", "Posterior 1", "Posterior 2"), each=m+1))
pdcilludf2$lk <- ordered(pdcilludf2$lk)
pdcillu2 <- ggplot(pdcilludf2, aes(lk, plk, group=which, linetype=which)) + geom_point() + geom_line() +
  bottomlegend + xlab(expression(l[k])) + ylab(expression(P(C[t]^k == l[k])))
pdcillu2

jointdf <- rbind(data.frame(bddf, which="Beta densities"),
                 data.frame(pdcilludf, which="Beta-binomial pmfs"))
joint <- ggplot(jointdf, aes(x, value, group=variable, linetype=variable)) + geom_point() + geom_line() +
#  facet_grid(. ~ which, scales="free") + bottomlegend
  facet_wrap(~which, scales="free") + rightlegend
joint # todo: axis labels, for beta densities only line not points

# better with grid.arrange:
# margins t, r, b, l
betadens2 <- betadens2 + theme(plot.margin = unit(c(0,0.5,0,0), "lines"))
pdcillu   <- pdcillu   + theme(plot.margin = unit(c(0,0,0,0), "lines"),
                               legend.margin = unit(-0.5, "cm"))
pdf("singleprior-pdc.pdf", width=8, height=3)
grid.arrange(betadens2, pdcillu, nrow=1, ncol=2, widths=c(1,1.2))
dev.off()

betapdfcdf      <- betapdfcdf    + theme(plot.margin = unit(c(0,0.5,0,-0.75), "lines"))
pdcillupmfcmf   <- pdcillupmfcmf + theme(plot.margin = unit(c(0,0,  0,-0.75), "lines"),
                                    legend.margin = unit(-0.5, "cm"))
pdf("singleprior-pdc2.pdf", width=8, height=5)
grid.arrange(betapdfcdf, pdcillupmfcmf, nrow=1, ncol=2, widths=c(1,1.2))
dev.off()


# ------------------- sets of Beta pdfs, parameter sets, sets of Beta-Binomial cmfs -------------------------

setpr <- BinomialLuckModel(n0=c(1,8), y0=c(0.7, 0.8))
setpos1 <- setpr
data(setpos1) <- BinomialData(s=12, n=16)
setpos2 <- setpr
data(setpos2) <- BinomialData(s= 0, n=16)

# Beta cdfs
betavec <- seq(0,1, length.out=201)

#cdfplot(setpos1, xvec=betavec)
#cdfplot(setpos1, xvec=betavec, control=controlList(posterior=TRUE))
#cdfplot(setpos2, xvec=betavec, control=controlList(posterior=TRUE))

betapriolow <- sapply(betavec, function(x){
  optim(par=c(4, 0.7), fn=function(.n0y0, x){
    pbetany(x, .n0y0)
  }, method="L-BFGS-B", control=list(fnscale=1),
  lower=c(n0(setpr)[1], y0(setpr)[1]),
  upper=c(n0(setpr)[2], y0(setpr)[2]), x=x)$value
})

betaprioupp <- sapply(betavec, function(x){
  optim(par=c(4, 0.8), fn=function(.n0y0, x){
    pbetany(x, .n0y0)
  }, method="L-BFGS-B", control=list(fnscale=-1),
  lower=c(n0(setpr)[1], y0(setpr)[1]),
  upper=c(n0(setpr)[2], y0(setpr)[2]), x=x)$value
})

betapriocorners <- sapply(luck4cny(setpr), function(x){
  pbetany(betavec, x)
})

betapos1low <- sapply(betavec, function(x){
  optim(par=c(4, 0.7), fn=function(.n0y0, x){
    pbetany(x, nyupdate(pr=.n0y0, data=c(tau(data(setpos1)), n(data(setpos1)))))
  }, method="L-BFGS-B", control=list(fnscale=1),
  lower=c(n0(setpr)[1], y0(setpr)[1]),
  upper=c(n0(setpr)[2], y0(setpr)[2]), x=x)$value
})

betapos1upp <- sapply(betavec, function(x){
  optim(par=c(4, 0.8), fn=function(.n0y0, x){
    pbetany(x, nyupdate(pr=.n0y0, data=c(tau(data(setpos1)), n(data(setpos1)))))
  }, method="L-BFGS-B", control=list(fnscale=-1),
  lower=c(n0(setpr)[1], y0(setpr)[1]),
  upper=c(n0(setpr)[2], y0(setpr)[2]), x=x)$value
})

betapos1corners <- sapply(luck4cny(setpos1, posterior=TRUE), function(x){
  pbetany(betavec, x)
})

betapos2low <- sapply(betavec, function(x){
  optim(par=c(4, 0.7), fn=function(.n0y0, x){
    pbetany(x, nyupdate(pr=.n0y0, data=c(tau(data(setpos2)), n(data(setpos2)))))
  }, method="L-BFGS-B", control=list(fnscale=1),
  lower=c(n0(setpr)[1], y0(setpr)[1]),
  upper=c(n0(setpr)[2], y0(setpr)[2]), x=x)$value
})

betapos2upp <- sapply(betavec, function(x){
  optim(par=c(4, 0.8), fn=function(.n0y0, x){
    pbetany(x, nyupdate(pr=.n0y0, data=c(tau(data(setpos2)), n(data(setpos2)))))
  }, method="L-BFGS-B", control=list(fnscale=-1),
  lower=c(n0(setpr)[1], y0(setpr)[1]),
  upper=c(n0(setpr)[2], y0(setpr)[2]), x=x)$value
})

betapos2corners <- sapply(luck4cny(setpos2, posterior=TRUE), function(x){
  pbetany(betavec, x)
})

betasetdf <- rbind(data.frame(x=betavec, Lower=betapriolow, Upper=betaprioupp, betapriocorners, Item="Prior", Facet="Prior & Posterior 1"),
                   data.frame(x=betavec, Lower=betapos1low, Upper=betapos1upp, betapos1corners, Item="Posterior 1", Facet="Prior & Posterior 1"),
                   data.frame(x=betavec, Lower=betapriolow, Upper=betaprioupp, betapriocorners, Item="Prior", Facet="Prior & Posterior 2"),
                   data.frame(x=betavec, Lower=betapos2low, Upper=betapos2upp, betapos2corners, Item="Posterior 2", Facet="Prior & Posterior 2"))
betasetdf$Item <- ordered(betasetdf$Item, levels=c("Prior", "Posterior 1", "Posterior 2"))

betaset1 <- ggplot(betasetdf) + geom_ribbon(aes(x=x, ymin=Lower, ymax=Upper, group=Item, colour=Item, fill=Item), alpha=0.3)
betaset1 <- betaset1 +
  geom_line(aes(x=x, y=c1, group=Item, colour=Item)) +
  geom_line(aes(x=x, y=c2, group=Item, colour=Item)) +
  geom_line(aes(x=x, y=c3, group=Item, colour=Item)) +
  geom_line(aes(x=x, y=c4, group=Item, colour=Item))
betaset1 <- betaset1 + facet_grid(Facet ~ .) + rightlegend + xlab(expression(p[t]^k)) + ylab("cdf")
betaset1

# parameter sets

n0vec <- seq(n0(setpr)[1], n0(setpr)[2], length.out=201)
prioparam <- data.frame(n0=n0vec, Lower=y0(setpr)[1], Upper=y0(setpr)[2])
pos1param <- data.frame(n0=updateLuckN(n0vec, n(data(setpos1))),
                        Lower=updateLuckY(n0vec, y0(setpos1)[1], tau(data(setpos1)), n(data(setpos1))),
                        Upper=updateLuckY(n0vec, y0(setpos1)[2], tau(data(setpos1)), n(data(setpos1))))
pos2param <- data.frame(n0=updateLuckN(n0vec, n(data(setpos2))),
                        Lower=updateLuckY(n0vec, y0(setpos2)[1], tau(data(setpos2)), n(data(setpos2))),
                        Upper=updateLuckY(n0vec, y0(setpos2)[2], tau(data(setpos2)), n(data(setpos2))))

paramsetdf <- rbind(data.frame(prioparam, Item="Prior", Facet="Prior & Posterior 1"),
                    data.frame(pos1param, Item="Posterior 1", Facet="Prior & Posterior 1"),
                    data.frame(prioparam, Item="Prior", Facet="Prior & Posterior 2"),
                    data.frame(pos2param, Item="Posterior 2", Facet="Prior & Posterior 2"))

paramset1 <- ggplot(paramsetdf) + geom_ribbon(aes(x=n0, ymin=Lower, ymax=Upper, group=Item, colour=Item, fill=Item), alpha=0.3)
paramset1 <- paramset1 + facet_grid(Facet ~ .) + guides(linetype="none", colour="none", fill="none")
paramset1 <- paramset1 + xlab(expression(n^(0))) + ylab(expression(y^(0)))
paramset1

paramset1 <- paramset1 + theme(plot.margin = unit(c(0,0.5,0,-0.25), "lines"))
betaset1  <- betaset1  + theme(plot.margin = unit(c(0,0,  0,-0.25), "lines"),
                                         legend.margin = unit(-0.5, "cm"))
# parameter sets and cdf sets in one plot
pdf("priorset.pdf", width=8, height=5)
grid.arrange(paramset1, betaset1, nrow=1, ncol=2, widths=c(1,1.2))
dev.off()

paramset2 <- ggplot(paramsetdf) + geom_ribbon(aes(x=n0, ymin=Lower, ymax=Upper, group=Item, colour=Item, fill=Item), alpha=0.3)
paramset2 <- paramset2 + xlab(expression(n^(0))) + ylab(expression(y^(0)))
paramset2

# plot with parameter sets only
pdf("paramsets.pdf", width=6, height=5)
paramset2 + theme(plot.margin = unit(c(0,0.5,0,-0.25), "lines"), legend.margin = unit(-0.5, "cm"))
dev.off()

# sets of Beta-Binomial cmfs
m <- 5
lvec <- 0:m

bebinpriolow <- sapply(lvec, function(x){
  optim(par=c(4, 0.7), fn=function(.n0y0, m, x){
    pbebinvec(.n0y0, m=m, l=x)
  }, method="L-BFGS-B", control=list(fnscale=1),
  lower=c(n0(setpr)[1], y0(setpr)[1]),
  upper=c(n0(setpr)[2], y0(setpr)[2]), m=m, x=x)$value
})

bebinprioupp <- sapply(lvec, function(x){
  optim(par=c(4, 0.8), fn=function(.n0y0, m, x){
    pbebinvec(.n0y0, m=m, l=x)
  }, method="L-BFGS-B", control=list(fnscale=-1),
  lower=c(n0(setpr)[1], y0(setpr)[1]),
  upper=c(n0(setpr)[2], y0(setpr)[2]), m=m, x=x)$value
})

bebinpriocorners <- sapply(luck4cny(setpr), function(x){
  pbebinvec(x, m=m)
})

bebinpos1low <- sapply(lvec, function(x){
  optim(par=c(4, 0.7), fn=function(.n0y0, m, x){
    pbebinvec(nyupdate(pr=.n0y0, data=c(tau(data(setpos1)), n(data(setpos1)))), m=m, l=x)
  }, method="L-BFGS-B", control=list(fnscale=1),
  lower=c(n0(setpr)[1], y0(setpr)[1]),
  upper=c(n0(setpr)[2], y0(setpr)[2]), m=m, x=x)$value
})

bebinpos1upp <- sapply(lvec, function(x){
  optim(par=c(4, 0.8), fn=function(.n0y0, m, x){
    pbebinvec(nyupdate(pr=.n0y0, data=c(tau(data(setpos1)), n(data(setpos1)))), m=m, l=x)
  }, method="L-BFGS-B", control=list(fnscale=-1),
  lower=c(n0(setpr)[1], y0(setpr)[1]),
  upper=c(n0(setpr)[2], y0(setpr)[2]), m=m, x=x)$value
})

bebinpos1corners <- sapply(luck4cny(setpr), function(x){
  pbebinvec(nyupdate(pr=x, data=c(tau(data(setpos1)), n(data(setpos1)))), m=m)
})

bebinpos2low <- sapply(lvec, function(x){
  optim(par=c(4, 0.7), fn=function(.n0y0, m, x){
    pbebinvec(nyupdate(pr=.n0y0, data=c(tau(data(setpos2)), n(data(setpos2)))), m=m, l=x)
  }, method="L-BFGS-B", control=list(fnscale=1),
  lower=c(n0(setpr)[1], y0(setpr)[1]),
  upper=c(n0(setpr)[2], y0(setpr)[2]), m=m, x=x)$value
})

bebinpos2upp <- sapply(lvec, function(x){
  optim(par=c(4, 0.8), fn=function(.n0y0, m, x){
    pbebinvec(nyupdate(pr=.n0y0, data=c(tau(data(setpos2)), n(data(setpos2)))), m=m, l=x)
  }, method="L-BFGS-B", control=list(fnscale=-1),
  lower=c(n0(setpr)[1], y0(setpr)[1]),
  upper=c(n0(setpr)[2], y0(setpr)[2]), m=m, x=x)$value
})

bebinpos2corners <- sapply(luck4cny(setpr), function(x){
  pbebinvec(nyupdate(pr=x, data=c(tau(data(setpos2)), n(data(setpos2)))), m=m)
})

bebinsetdf <- rbind(data.frame(x=lvec, Lower=bebinpriolow, Upper=bebinprioupp, bebinpriocorners, Item="Prior", Facet="Prior & Posterior 1"),
                    data.frame(x=lvec, Lower=bebinpos1low, Upper=bebinpos1upp, bebinpos1corners, Item="Posterior 1", Facet="Prior & Posterior 1"),
                    data.frame(x=lvec, Lower=bebinpriolow, Upper=bebinprioupp, bebinpriocorners, Item="Prior", Facet="Prior & Posterior 2"),
                    data.frame(x=lvec, Lower=bebinpos2low, Upper=bebinpos2upp, bebinpos2corners, Item="Posterior 2", Facet="Prior & Posterior 2"))
bebinsetdf$Item <- ordered(bebinsetdf$Item, levels=c("Prior", "Posterior 1", "Posterior 2"))

#
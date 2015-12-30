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

# ------------------- sets of Beta pdfs, sets of Beta-Bin cmfs -------------------------

setpr <- BinomialLuckModel(n0=c(1,8), y0=c(0.7, 0.8))
setpos1 <- setpr
data(setpos1) <- BinomialData(s=12, n=16)
setpos2 <- setpr
data(setpos2) <- BinomialData(s= 0, n=16)

betavec <- seq(0,1, length.out=200)

cdfplot(setpos1, xvec=betavec)
cdfplot(setpos1, xvec=betavec, control=controlList(posterior=TRUE))
cdfplot(setpos2, xvec=betavec, control=controlList(posterior=TRUE))

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

betasetdf <- rbind(data.frame(x=betavec, Lower=betapriolow, Upper=betaprioupp, Item="Prior"),
                   data.frame(x=betavec, Lower=betapos1low, Upper=betapos1upp, Item="Posterior 1"),
                   data.frame(x=betavec, Lower=betapos2low, Upper=betapos2upp, Item="Posterior 2"))
betasetdf$Item <- ordered(betasetdf$Item, levels=c("Prior", "Posterior 1", "Posterior 2"))
betaset1 <- ggplot(betasetdf) + geom_ribbon(aes(x=x, ymin=Lower, ymax=Upper, group=Item, colour=Item, fill=Item), alpha=0.3)
betaset1

#
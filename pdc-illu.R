# Illustration prior-data conflict in Beta-Binomial model
# plots of posterior predictive

library(ggplot2)
library(reshape2)
library(grid)
library(gridExtra)

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

# beta densities
betavec <- seq(0,1, length.out=200)
bdsgpr <- dbetany(betavec, sgpr)
bdpos1 <- dbetany(betavec, pos1)
bdpos2 <- dbetany(betavec, pos2)

bottomlegend <- theme(legend.position = 'bottom', legend.direction = 'horizontal', legend.title = element_blank())
rightlegend <- theme(legend.title = element_blank())

bddf <- data.frame(x=betavec, Prior=bdsgpr, "Posterior1"=bdpos1, "Posterior2"=bdpos2)
names(bddf)[c(3,4)] <- c("Posterior 1", "Posterior 2")
bddf <- melt(bddf, "x")

betadens <- ggplot(bddf, aes(x, value, group=variable, linetype=variable)) + geom_line() +
  bottomlegend + xlab(expression(p[t]^k)) + ylab(expression(f(p[t]^k)))
betadens

betadens2 <- ggplot(bddf, aes(x, value, group=variable, linetype=variable)) + geom_line() +
  guides(linetype="none") + xlab(expression(p[t]^k)) + ylab(expression(f(p[t]^k)))
betadens2

# beta-binomial distributions
#dbebinvec(pos1,2,1)
m <- 5
bebinpr <- dbebinvec(sgpr, m)
bebin1 <- dbebinvec(pos1, m)
bebin2 <- dbebinvec(pos2, m)

pdcilludf <- data.frame(x=rep(0:m), Prior=bebinpr, "Posterior1"=bebin1, "Posterior2"=bebin2)
names(pdcilludf)[c(3,4)] <- c("Posterior 1", "Posterior 2")
pdcilludf$x <- ordered(pdcilludf$x)
pdcilludf <- melt(pdcilludf, "x")
pdcillu <- ggplot(pdcilludf, aes(x, value, group=variable, linetype=variable)) + geom_point() + geom_line() +
  rightlegend + xlab(expression(l[k])) + ylab(expression(P(C[t]^k == l[k])))
pdcillu

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

pdf("test.pdf")
#betadens2
pdcillu
dev.off()

#
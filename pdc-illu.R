# Illustration prior-data conflict in Beta-Binomial model
# plots of posterior predictive

library(ggplot2)
library(reshape2)

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

bddf <- melt(data.frame(x=betavec, Prior=bdsgpr, "Posterior 1"=bdpos1, "Posterior 2"=bdpos2), "x")

betadens <- ggplot(bddf, aes(x, value, group=variable, linetype=variable)) + geom_line() +
  bottomlegend + xlab(expression(p[t]^k)) + ylab(expression(f(p[t]^k)))
betadens


# beta-binomial distributions
#dbebinvec(pos1,2,1)
m <- 4
bebinpr <- dbebinvec(sgpr, m)
bebin1 <- dbebinvec(pos1, m)
bebin2 <- dbebinvec(pos2, m)

bottomlegend <- theme(legend.position = 'bottom', legend.direction = 'horizontal', legend.title = element_blank())
rightlegend <- theme(legend.title = element_blank())

pdcilludf <- data.frame(lk=rep(0:m, 2), plk=c(bebin1, bebin2), which=rep(c("Posterior 1", "Posterior 2"), each=m+1))

pdcillu <- ggplot(pdcilludf, aes(lk, plk, group=which, linetype=which)) + geom_point() + geom_line() +
  bottomlegend + xlab(expression(l[k])) + ylab(expression(P(C[t]^k == l[k])))
pdcillu

pdcillu <- ggplot(pdcilludf, aes(lk, plk)) + geom_point() + geom_line() + facet_grid(. ~ which) +
  xlab(expression(l[k])) + ylab(expression(P(C[t]^k == l[k]))) + 
  scale_x_discrete(breaks=0:m) + coord_cartesian(xlim=c(-0.5,m+0.5)) # axis tics at 0:m
pdcillu

#
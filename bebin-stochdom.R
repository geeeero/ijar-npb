dbebinvec <- function(pos, m, l = NULL){
  if (is.null(l))
    l <- 0:m
  sapply(l, function(li, m, n, y){
    choose(m,li) * beta(li+n*y, m-li+n*(1-y)) / beta(n*y, n*(1-y)) 
  }, m=m, n=pos[1], y=pos[2])
}

pbebinvec <- function(...)
  cumsum(dbebinvec(...))

# beta-binomial is stochastically ordered along y:
yvec <- (1:9)/10
pbebinmat <- sapply(yvec, function(y,n,m){
  pbebinvec(c(n,y),m=m)
}, n=1, m=5)
matplot(0:5, pbebinmat, type="l")
all(apply(pbebinmat[-6,], 1, order) == matrix(9:1, nrow=9, ncol=5))

# ...but not along n:
nvec <- 1:9/10
pbebinmat <- sapply(nvec, function(n,y,m){
  pbebinvec(c(n,y),m=m)
}, y=0.5, m=5)
matplot(0:5, pbebinmat, type="l")

#
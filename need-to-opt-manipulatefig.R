p <- seq(0, 1, length.out = 100)

prob.s <- function(y, p, N, m) {
  pbinom(floor(y*(N+m-1)), N, p) - pbinom(ceiling(y*(N+m-1)-m+1) - 1, N, p)
}
prob <- Vectorize(prob.s, "p")

prob2.s <- function(y, p, N, m, n.u, n.l) {
  s.bad <- NULL
  for(s in 0:N) {
    L0 <- (prod(0:(m-1)+n.u*(1-y)+N-s)*prod(0:(m-1)+n.l+N)) /
      (prod(0:(m-1)+n.l*(1-y)+N-s)*prod(0:(m-1)+n.u+N))
    Lm <- (prod(0:(m-1)+n.u*y+s)*prod(0:(m-1)+n.l+N)) /
      (prod(0:(m-1)+n.l*y+s)*prod(0:(m-1)+n.u+N))
    if(is.nan(Lm) || is.nan(L0) || (!(L0<=1 && Lm>=1) && !(L0>=1 && Lm<=1))) {
      s.bad <- c(s.bad, s)
    }
  }
  if(length(s.bad)==0)
    return(0)
  sum(dbinom(s.bad, N, p))
}
prob2 <- Vectorize(prob2.s, "p")

manipulate(
{
  plot(p, prob(y, p, N, m), type="l", ylim=c(0,1), ylab="Probability must optimise")
  points(p, prob2(y, p, N, m, n.u, n.l), type="l", col="red")
  abline(v = y)
},
y = slider(0, 1, step=0.1, initial = 0.1),
N = slider(5, 100, step=5, initial = 20),
m = slider(1, 10, initial = 3),
n.l = slider(1, 20),
n.u = slider(1, 20, initial = 5)
)


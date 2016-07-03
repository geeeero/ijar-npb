library(ggplot2)

p <- seq(0, 1, length.out = 100)

prob.s <- function(y, p, N, m) {
  if(p<0)
    return(NA)
  if(p>1)
    return(NA)
  pbinom(floor(y*(N+m-1)), N, p) - pbinom(ceiling(y*(N+m-1)-m+1) - 1, N, p)
}
prob <- Vectorize(prob.s, "p")

prob2.s <- function(y, p, N, m, n.u, n.l) {
  if(p<0)
    return(NA)
  if(p>1)
    return(NA)

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

# Plot 1
y <- 0.1; N <- 100; m <- 3; n.l <- 1; n.u <- 5

y <- 0.1
d1 <- data.frame(p=rep(p, 2),
                 popt=c(prob(y, p, N, m), prob2(y, p, N, m, n.u, n.l)),
                 lab=rep(c("Theorem 2", "Lemma 3"), each=length(p)),
                 y=y)
y <- 0.25
d2 <- data.frame(p=rep(p, 2),
                 popt=c(prob(y, p, N, m), prob2(y, p, N, m, n.u, n.l)),
                 lab=rep(c("Theorem 2", "Lemma 3"), each=length(p)),
                 y=y)
y <- 0.5
d3 <- data.frame(p=rep(p, 2),
                 popt=c(prob(y, p, N, m), prob2(y, p, N, m, n.u, n.l)),
                 lab=rep(c("Theorem 2", "Lemma 3"), each=length(p)),
                 y=y)


pdf("need-to-opt-1.pdf", width=8, height=2)
ggplot(rbind(d1, d2, d3)) +
  ylim(0, 1) + xlab("p") + ylab("Prob need to optimise") + theme_bw() +
  scale_linetype_discrete(name="Theory", breaks=c("Theorem 2", "Lemma 3")) +
  facet_grid(.~y, labeller = label_both) +
  geom_line(aes(x=p, y=popt, lty=lab))
dev.off()

# Plot 2
y <- 0.1; N <- 20; m <- 3; n.l <- 1; n.u <- 5

N <- 10
d1 <- data.frame(p=rep(p, 2),
                 popt=c(prob(y, p, N, m), prob2(y, p, N, m, n.u, n.l)),
                 lab=rep(c("Theorem 2", "Lemma 3"), each=length(p)),
                 N=10)
N <- 100
d2 <- data.frame(p=rep(p, 2),
                 popt=c(prob(y, p, N, m), prob2(y, p, N, m, n.u, n.l)),
                 lab=rep(c("Theorem 2", "Lemma 3"), each=length(p)),
                 N=100)
N <- 1000
d3 <- data.frame(p=rep(p, 2),
                 popt=c(prob(y, p, N, m), prob2(y, p, N, m, n.u, n.l)),
                 lab=rep(c("Theorem 2", "Lemma 3"), each=length(p)),
                 N=1000)


pdf("need-to-opt-2.pdf", width=8, height=2)
ggplot(rbind(d1, d2, d3)) +
  ylim(0, 1) + xlab("p") + ylab("Prob need to optimise") + theme_bw() +
  scale_linetype_discrete(name="Theory", breaks=c("Theorem 2", "Lemma 3")) +
  facet_grid(.~N, labeller = label_both) +
  geom_line(aes(x=p, y=popt, lty=lab))
dev.off()

# Plot 3
y <- 0.4; N <- 100; m <- 3; n.l <- 1; n.u <- 5

m <- 1
d1 <- data.frame(p=rep(p, 2),
                 popt=c(prob(y, p, N, m), prob2(y, p, N, m, n.u, n.l)),
                 lab=rep(c("Theorem 2", "Lemma 3"), each=length(p)),
                 m=1)
m <- 5
d2 <- data.frame(p=rep(p, 2),
                 popt=c(prob(y, p, N, m), prob2(y, p, N, m, n.u, n.l)),
                 lab=rep(c("Theorem 2", "Lemma 3"), each=length(p)),
                 m=5)
m <- 10
d3 <- data.frame(p=rep(p, 2),
                 popt=c(prob(y, p, N, m), prob2(y, p, N, m, n.u, n.l)),
                 lab=rep(c("Theorem 2", "Lemma 3"), each=length(p)),
                 m=10)


pdf("need-to-opt-3.pdf", width=8, height=2)
ggplot(rbind(d1, d2, d3)) +
  ylim(0, 1) + xlab("p") + ylab("Prob need to optimise") + theme_bw() +
  scale_linetype_discrete(name="Theory", breaks=c("Theorem 2", "Lemma 3")) +
  facet_grid(.~m, labeller = label_both) +
  geom_line(aes(x=p, y=popt, lty=lab))
dev.off()

# Extra function defs
prob.p <- Vectorize(prob.s, c("y","p"))
prob2.p <- Vectorize(prob2.s, c("y", "p"))

# Plot 4
y <- seq(0, 1, length.out = 51)

N <- 100; m <- 3; n.l <- 1; n.u <- 5

N <- 50
d1 <- data.frame(p=rep(y, 2),
                 popt=c(prob.p(y, y, N, m), prob2.p(y, y, N, m, n.u, n.l)),
                 lab=rep(c("Theorem 2", "Lemma 3"), each=length(y)),
                 N="10")
N <- 100
d2 <- data.frame(p=rep(y, 2),
                 popt=c(prob.p(y, y, N, m), prob2.p(y, y, N, m, n.u, n.l)),
                 lab=rep(c("Theorem 2", "Lemma 3"), each=length(y)),
                 N="100")
N <- 1000
d3 <- data.frame(p=rep(y, 2),
                 popt=c(prob.p(y, y, N, m), prob2.p(y, y, N, m, n.u, n.l)),
                 lab=rep(c("Theorem 2", "Lemma 3"), each=length(y)),
                 N="1000")


pdf("need-to-opt-4.pdf", width=3.95, height=3)
ggplot() +
  ylim(0, 1) + xlab("y") + ylab("Worst case prob need to optimise") + theme_bw() +
  scale_linetype_discrete(name="Theory", breaks=c("Theorem 2", "Lemma 3")) +
  geom_line(data=d1, aes(x=p, y=popt, lty=lab, colour=N)) +
  geom_line(data=d2, aes(x=p, y=popt, lty=lab, colour=N)) +
  geom_line(data=d3, aes(x=p, y=popt, lty=lab, colour=N))
dev.off()

# Plot 5
y <- seq(0, 1, length.out = 51)

N <- 100; m <- 3; n.l <- 1; n.u <- 5

m <- 2
d1 <- data.frame(p=rep(y, 2),
                 popt=c(prob.p(y, y, N, m), prob2.p(y, y, N, m, n.u, n.l)),
                 lab=rep(c("Theorem 2", "Lemma 3"), each=length(y)),
                 m="2")
m <- 5
d2 <- data.frame(p=rep(y, 2),
                 popt=c(prob.p(y, y, N, m), prob2.p(y, y, N, m, n.u, n.l)),
                 lab=rep(c("Theorem 2", "Lemma 3"), each=length(y)),
                 m="5")
m <- 10
d3 <- data.frame(p=rep(y, 2),
                 popt=c(prob.p(y, y, N, m), prob2.p(y, y, N, m, n.u, n.l)),
                 lab=rep(c("Theorem 2", "Lemma 3"), each=length(y)),
                 m="10")


pdf("need-to-opt-5.pdf", width=3.95, height=3)
ggplot() +
  ylim(0, 1) + xlab("y") + ylab("Worst case prob need to optimise") + theme_bw() +
  scale_linetype_discrete(name="Theory", breaks=c("Theorem 2", "Lemma 3")) +
  scale_colour_discrete(breaks=c("2", "5", "10")) +
  geom_line(data=d1, aes(x=p, y=popt, lty=lab, colour=m)) +
  geom_line(data=d2, aes(x=p, y=popt, lty=lab, colour=m)) +
  geom_line(data=d3, aes(x=p, y=popt, lty=lab, colour=m))
dev.off()

# Plot 6
y <- seq(0, 1, length.out = 51)

N <- 100; m <- 3; n.l <- 1; n.u <- 5

N <- 50
d1 <- data.frame(p=rep(y, 2),
                 popt=c(pmax(prob.p(y, y-0.1, N, m), prob.p(y, y+0.1, N, m), na.rm = TRUE),
                        pmax(prob2.p(y, y-0.1, N, m, n.u, n.l), prob2.p(y, y+0.1, N, m, n.u, n.l), na.rm = TRUE)),
                 lab=rep(c("Theorem 2", "Lemma 3"), each=length(y)),
                 N="10")
N <- 100
d2 <- data.frame(p=rep(y, 2),
                 popt=c(pmax(prob.p(y, y-0.1, N, m), prob.p(y, y+0.1, N, m), na.rm = TRUE),
                        pmax(prob2.p(y, y-0.1, N, m, n.u, n.l), prob2.p(y, y+0.1, N, m, n.u, n.l), na.rm = TRUE)),
                 lab=rep(c("Theorem 2", "Lemma 3"), each=length(y)),
                 N="100")
N <- 1000
d3 <- data.frame(p=rep(y, 2),
                 popt=c(pmax(prob.p(y, y-0.1, N, m), prob.p(y, y+0.1, N, m), na.rm = TRUE),
                        pmax(prob2.p(y, y-0.1, N, m, n.u, n.l), prob2.p(y, y+0.1, N, m, n.u, n.l), na.rm = TRUE)),
                 lab=rep(c("Theorem 2", "Lemma 3"), each=length(y)),
                 N="1000")


pdf("need-to-opt-6.pdf", width=3.95, height=3)
ggplot() +
  ylim(0, 1) + xlab("y") + ylab("Worst case prob need to optimise") + theme_bw() +
  scale_linetype_discrete(name="Theory", breaks=c("Theorem 2", "Lemma 3")) +
  geom_line(data=d1, aes(x=p, y=popt, lty=lab, colour=N)) +
  geom_line(data=d2, aes(x=p, y=popt, lty=lab, colour=N)) +
  geom_line(data=d3, aes(x=p, y=popt, lty=lab, colour=N))
dev.off()

# Plot 7
y <- seq(0, 1, length.out = 51)

N <- 100; m <- 3; n.l <- 1; n.u <- 5

m <- 2
d1 <- data.frame(p=rep(y, 2),
                 popt=c(pmax(prob.p(y, y-0.1, N, m), prob.p(y, y+0.1, N, m), na.rm = TRUE),
                        pmax(prob2.p(y, y-0.1, N, m, n.u, n.l), prob2.p(y, y+0.1, N, m, n.u, n.l), na.rm = TRUE)),
                 lab=rep(c("Theorem 2", "Lemma 3"), each=length(y)),
                 m="2")
m <- 5
d2 <- data.frame(p=rep(y, 2),
                 popt=c(pmax(prob.p(y, y-0.1, N, m), prob.p(y, y+0.1, N, m), na.rm = TRUE),
                        pmax(prob2.p(y, y-0.1, N, m, n.u, n.l), prob2.p(y, y+0.1, N, m, n.u, n.l), na.rm = TRUE)),
                 lab=rep(c("Theorem 2", "Lemma 3"), each=length(y)),
                 m="5")
m <- 10
d3 <- data.frame(p=rep(y, 2),
                 popt=c(pmax(prob.p(y, y-0.1, N, m), prob.p(y, y+0.1, N, m), na.rm = TRUE),
                        pmax(prob2.p(y, y-0.1, N, m, n.u, n.l), prob2.p(y, y+0.1, N, m, n.u, n.l), na.rm = TRUE)),
                 lab=rep(c("Theorem 2", "Lemma 3"), each=length(y)),
                 m="10")


pdf("need-to-opt-7.pdf", width=3.95, height=3)
ggplot() +
  ylim(0, 1) + xlab("y") + ylab("Worst case prob need to optimise") + theme_bw() +
  scale_linetype_discrete(name="Theory", breaks=c("Theorem 2", "Lemma 3")) +
  scale_colour_discrete(breaks=c("2", "5", "10")) +
  geom_line(data=d1, aes(x=p, y=popt, lty=lab, colour=m)) +
  geom_line(data=d2, aes(x=p, y=popt, lty=lab, colour=m)) +
  geom_line(data=d3, aes(x=p, y=popt, lty=lab, colour=m))
dev.off()

#### Brake example ####
t <- seq(0, 10, length.out=301)

# Component M
priorU <- 1-pweibull(t, shape=2.5, scale=8)
priorU[priorU==1] <- 1-1e-6
priorL <- 1-pweibull(t, shape=2.5, scale=6)
priorL[priorL==1] <- 1-1e-6

m <- 1; N <- 5; nl <- 1; nu <- 8

M <- data.frame(t=rep(t, 2),
                popt=c(pmax(prob2.p(priorU, priorU-0.1, N, m, n.u, n.l), prob2.p(priorU, priorU+0.1, N, m, n.u, n.l), na.rm = TRUE),
                       pmax(prob2.p(priorL, priorL-0.1, N, m, n.u, n.l), prob2.p(priorL, priorL+0.1, N, m, n.u, n.l), na.rm = TRUE)),
                lab=rep(c("Upper", "Lower"), each=length(t)),
                comp="M")

# Component H
priorU <- rep(0.9999, 301)
priorL <- rep(0.0001, 301)

m <- 1; N <- 10; nl <- 1; nu <- 2

H <- data.frame(t=rep(t, 2),
                popt=c(pmax(prob2.p(priorU, priorU-0.1, N, m, n.u, n.l), prob2.p(priorU, priorU+0.1, N, m, n.u, n.l), na.rm = TRUE),
                       pmax(prob2.p(priorL, priorL-0.1, N, m, n.u, n.l), prob2.p(priorL, priorL+0.1, N, m, n.u, n.l), na.rm = TRUE)),
                lab=rep(c("Upper", "Lower"), each=length(t)),
                comp="H")

# Component C
priorU <- c(rep(c(0.99, 0.98, 0.96, 0.95, 0.90, 0.65, 0.50, 0.45, 0.43, 0.42), each=30), 0.42)
priorL <- c(rep(c(0.75, 0.73, 0.71, 0.70, 0.60, 0.45, 0.30, 0.23, 0.21, 0.20), each=30), 0.20)

m <- 4; N <- 15; nl <- 1; nu <- 2

C <- data.frame(t=rep(t, 2),
                popt=c(pmax(prob2.p(priorU, priorU-0.1, N, m, n.u, n.l), prob2.p(priorU, priorU+0.1, N, m, n.u, n.l), na.rm = TRUE),
                       pmax(prob2.p(priorL, priorL-0.1, N, m, n.u, n.l), prob2.p(priorL, priorL+0.1, N, m, n.u, n.l), na.rm = TRUE)),
                lab=rep(c("Upper", "Lower"), each=length(t)),
                comp="C")

# Component P
priorU <- c(rep(c(0.99, 0.65), each=150), 0.65)
priorL <- c(rep(c(0.5, 0.01), each=150), 0.01)

m <- 4; N <- 20; nl <- 1; nu <- 2

P <- data.frame(t=rep(t, 2),
                popt=c(pmax(prob2.p(priorU, priorU-0.1, N, m, n.u, n.l), prob2.p(priorU, priorU+0.1, N, m, n.u, n.l), na.rm = TRUE),
                       pmax(prob2.p(priorL, priorL-0.1, N, m, n.u, n.l), prob2.p(priorL, priorL+0.1, N, m, n.u, n.l), na.rm = TRUE)),
                lab=rep(c("Upper", "Lower"), each=length(t)),
                comp="P")

pdf("need-to-opt-8.pdf", width=8, height=4)
ggplot() +
  ylim(0, 1) + xlab("Time") + ylab("Worst case prob need to optimise") + theme_bw() +
  scale_linetype_discrete(name="Prior", breaks=c("Upper", "Lower")) +
  facet_wrap(~comp) +
  geom_line(data=M, aes(x=t, y=popt, lty=lab)) +
  geom_line(data=H, aes(x=t, y=popt, lty=lab)) +
  geom_line(data=C, aes(x=t, y=popt, lty=lab)) +
  geom_line(data=P, aes(x=t, y=popt, lty=lab))
dev.off()

library(dplyr)
mean(c(mean(t(M %>% filter(lab=="Lower") %>% select(popt))),
       mean(t(H %>% filter(lab=="Lower") %>% select(popt))),
       mean(t(C %>% filter(lab=="Lower") %>% select(popt))),
       mean(t(P %>% filter(lab=="Lower") %>% select(popt)))))

mean(c(mean(t(M %>% filter(lab=="Upper") %>% select(popt))),
       mean(t(H %>% filter(lab=="Upper") %>% select(popt))),
       mean(t(C %>% filter(lab=="Upper") %>% select(popt))),
       mean(t(P %>% filter(lab=="Upper") %>% select(popt)))))

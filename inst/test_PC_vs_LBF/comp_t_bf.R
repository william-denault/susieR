# Short script to compare different Bayes factor calculations.

# Compute Wakefield's asymptotic (log) Bayes factor.
compute_log_abf <- function (bhat, shat, s0)
  dnorm(bhat,0,sqrt(shat^2 + s0^2),log = TRUE) -
  dnorm(bhat,0,sqrt(shat^2),log = TRUE)

# This is William's code.
compute_log_ssbf_william <- function (x, y, s0, sd = 1) {
  n        <- length(x)
  X        <- cbind(1,x,1)
  invnu    <- diag(c(0,1/s0^2,1/sd^2))
  invOmega <- invnu + t(X) %*% X
  B        <- solve(invOmega,t(X) %*% cbind(y))
  return(-log(det(invOmega))/2 + log(n)/2 - log(s0)
         -log(sd) - n/2*(log(t(y - X %*% B) %*% y) -
                           log(t(y) %*% y - n*mean(y)^2)))
}

# This is the Bayes factor in which we treat the residual variance as
# unknown with a IG(a/2,b/2) prior, then set a = 0,b = 0. Here, s0^2
# is the prior variance of the coefficient, b (which is scaled by the
# residual variance).
compute_log_ssbf <- function (x, y, s0) {
  x   <- x - mean(x)
  y   <- y - mean(y)
  n   <- length(x)
  xx  <- sum(x*x)
  xy  <- sum(x*y)
  yy  <- sum(y*y)
  r0  <- s0/(s0 + 1/xx)
  sxy <- xy/sqrt(xx*yy)
  return((log(1 - r0) - n*log(1 - r0*sxy^2))/2)
}

# Compute William's "ratio of t's" (log) Bayes factor.
compute_log_tbf <- function (bhat,shat,s0,df)
  dst(bhat,0,sqrt(shat^2 + s0^2),df,log = TRUE) -
  dst(bhat,0,shat,df,log = TRUE)

library(LaplacesDemon)
set.seed(1)

s0     <- 10
nr     <- 10000
w_lbf  <- rep(0,nr)
x_lbf  <- rep(0,nr)
ss_lbf <- rep(0,nr)
t_lbf  <- rep(0,nr)
for (i in 1:nr) {

  # Simulate a data set.
  n <- sample(c(10,20,30),1)
  b <- sample(c(-1,-0.1,0,0.1,1),1)
  x <- rnorm(n)
  x <- as.vector(scale(x,center = TRUE,scale = TRUE))
  y <- b*x + rnorm(n)
  y <- as.vector(scale(y,center = TRUE,scale = TRUE))

  # Compute the (log) Bayes factors.
  fit       <- lm(y ~ x - 1)
  b_mle     <- summary(fit)$coefficients
  # > sqrt(mean((y - x*b_mle["x","Estimate"])^2)*n/(n-1)) - summary(fit)$sigma
  # [1] -1.11e-16
  bhat      <- b_mle["x","Estimate"]
  shat      <- b_mle["x","Std. Error"]
  w_lbf[i]  <- compute_log_abf(bhat,shat,s0)
  x_lbf[i]  <- compute_log_ssbf_william(x,y,s0)
  ss_lbf[i] <- compute_log_ssbf(x,y,s0)
  t_lbf[i]  <- compute_log_tbf(bhat,shat,s0,n)
}

print(range(x_lbf - ss_lbf))
plot(ss_lbf,w_lbf,pch = 20,cex = 0.8,col = "darkblue",
     xlab = "exact logBF",ylab = "ratio-of-t logBF")
points(ss_lbf,t_lbf,pch = 20,cex = 0.8,col = "dodgerblue")
abline(a = 0,b = 1,col = "magenta",lty = "dotted")
print(cor(ss_lbf,t_lbf))

# Offset in approximation (when s0 is large).
print(coef(lm(ss_lbf ~ t_lbf)))
print(coef(lm( ss_lbf ~x_lbf)))

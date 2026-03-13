GaussApprox <- function(f,
                        n,
                        Q,
                        k = 5,
                        h = 0.005,
                        x = NULL) {
    dsmall <- .Machine$double.eps^0.5
    m <- ncol(Q)
    if (is.null(x))
        x <- double(m)
    for (j in 1:k) {
        fx <- f(x[1:n])
        fa  <- f(x[1:n] - h)
        fb  <- f(x[1:n] + h)
        cc <- (2*fx -fa -fb)/(h^2)
        cc[cc<dsmall] <- dsmall
        bb <- (fb - fa)/(2*h) + x[1:n]*cc
        Qn <- Q + Diagonal(m, c(cc, rep(0, m-n)))
        L <- chol(Qn) ### "aparentlty" L is not use, but IS
        x <- drop(solve(Qn, c(bb, rep(0, m-n))))
    }
    return(list(x=x, bb=bb, cc=cc, Q=Qn))
}

n <- 70
x <- sin(2*pi*1:n/n)

plot(x)

 y <- rpois(n, exp(1+x))
 par(mfrow=c(1,1), mar=c(2,3,0.5,0.5), mgp=c(2,0.5,0))
 plot(y, pch=19)
 ## define the likeihood function
  llp <- function(e)
   dpois(y, exp(e), log=TRUE)
## define the precision matrix structure

library(Matrix)
Q2 <- crossprod(diff(Diagonal(n), differences=2))

Q2[1:10, 1:10]

Qb <- Q2[1:10, 1:10]  + Diagonal(10)
Qb

str(Qb)

L <- chol(Qb)

str(Qb)

ga <- GaussApprox(llp, n, Q2*5000)

names(ga)

ga$V <- solve(ga$Q)
ga$sd <- sqrt(diag(ga$V))

tail(ga$x,2)

plot(x, type = 'l', lwd = 2)
lines(ga$x, col = 2, lty = 2)

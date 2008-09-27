`BreakMeans` <-
function(x){
# Find breaks and thereafter the means
#
# x  : series
#
    dwt.x <- modwt(x, wf = "haar", n.levels = 2)
    w1 <- as.numeric(dwt.x[["d1"]])
    w1<-round(w1,8) # rounding
    n <- length(w1)
    sig <- mad(w1)
    # here the constant is  0.6744908
    cn <- sig * (2 * log(n))^0.5
    i <- 1:n
    a <- i[abs(w1) > cn]
    if(length(a) != 0) {
       n.a <- length(a)
       g <- rep(0, n.a)
       k <- rep(0, n.a)
       s1 <- rep(0, n.a)
       for(j in 1:n.a) {
          s1[j] <- sign(w1[a[j]])
          y <- x[ - a[j]]
          # now do wavelet analysis and find new w1
          dwt.y <- modwt(y, wf = "haar", n.levels = 2)
          w1n <- as.numeric(dwt.y[["d1"]])
          w1n<-round(w1n,8)
          nn <- length(w1n)
          sig.n <- mad(w1n)
          # here the constant is  0.6744908
          cnn = sig.n * (2 * log(nn))^0.5
          # If any one of the following gives T, then the signal, that max(abs(w[j])) gives, is
          # a real signal of breakpoint.
          g[j] <- any((abs(w1n[(a[j]):(a[j])]) > cnn) & (sign(w1n[(a[j]):(a[j])]) == s1[j]))
          # forward
          k[j] <- any((abs(w1n[(a[j] - 1):(a[j] - 1)]) > cnn) & (sign(w1n[(a[j] - 1):(a[j] - 1)]) == s1[j]))
            }
        b.n <- (g == 1) | (k == 1)
        # It tests whether this is a break point or not. If T, then breakpoint, otherwise
        # not a break point
        br.pt <- a[b.n]
        # these are the break points
        diff.pt <- c(0, diff(br.pt))
        diff.pt.n <- br.pt[diff.pt != 1]
     } else  diff.pt.n <- NA
    bp <- c(1, diff.pt.n, length(x))
    # Need to consider the first and last values
    bp.n <- Unique(c(abs(bp), abs(diff.pt.n)))
    # Find whether the first and last values
    # appeared before. 
    bp.n <- as.numeric(na.exclude(bp.n))
    # If no jump point detected, then the result would provide only "NA". 
    g.mean <- rep(NA, length(bp.n))
    g.mean<- MeansBreakPoints(x, bp.n)
    g.mean<- c(g.mean, g.mean[length(g.mean)])
    all.mean <- rep(g.mean, c(diff(bp.n), 1))
    all.mean
}


`MeansBreakPoints` <-
function(x, bp){
#
# Finding the region means given the series and the break points. 
#
# x : Given series
# bp: Break points
#
    m <- length(bp)
    BMean <- rep(NA, m-1)
    if (m==2)  BMean<- mean(x[bp[1]:bp[2]])
    if (m>2) {
        for ( i in 1:(m-2)) BMean[i] <- mean(x[bp[i]:(bp[i+1]-1)])
            BMean[m-1]<- mean(x[bp[m-1]:bp[m]])
    }
     BMean
 }


`SymmetricMA` <-
function(x, w=5){
#
# Finding moving average smoothing using neighboring size w. 
#
# x : the series
# w : window size
#
    n<- length(x)
    z<-rep(0, n)
    for (k in (w+1):(n-w)) 
        z[k]<- (1/(2*w+1))* sum(x[(k-w):(k+w)])

    for (k in 1:w){
        s<-max(1, k-w)
        t<-min(n, k+w)
        z[k]<- sum(x[s:t])/(t-s+1)
        }   

    for (k in (n-w+1):n){
        s<-max(1, k-w)
        t<-min(n, k+w)
        z[k]<- sum(x[s:t])/(t-s+1)
        }   
    z   
  }


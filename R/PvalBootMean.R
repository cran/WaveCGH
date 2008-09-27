`PvalBootMean` <-
function(x, ModelSelection=FALSE, nsim=200){
#
# Finds p-value for testing whether the mean is significantly different from zero.
#
# x              : CGH data
# ModelSelection : Whether to select a model or just use AR(1) process
# nsim           : number of bootstrap simulations 
    n<- length(x)
    Means<- BreakMeans(x)
    e<- x-Means
    if (ModelSelection){
        p <- ifelse(n>=8, SelectModel(e, lag.max=5, Criterion = "BIC", Best=1),1)
            } else p<-1
    est<- FitAR(e,p)
    sigA<- sqrt(est$sigsqHat)
    phiHat<- est$phiHat
    i<-1:n
    bp <- i[c(1, diff(Means)) != 0]
    bp <- Unique(c(bp, i[length(i)]))
    bp <- as.numeric(na.exclude(bp))

    e.boot<- apply(matrix(rep(1, nsim),ncol=1), 1, function(x)SimulateGaussianAR( phi=phiHat, n=n,InnovationVariance=sigA))

    OriginalMean<- MeansBreakPoints(x, bp)
    OriginalMeanMatrix<-matrix(rep(OriginalMean,nsim),ncol=nsim)
    BootMean <- apply(e.boot, 2, MeansBreakPoints, bp=bp)

    TFMatrix<- abs(BootMean)>abs(OriginalMeanMatrix)
    pval<-apply(TFMatrix, 1, function(x) sum(x)/nsim)
    pval.1<- c(pval, pval[length(pval)])
    AllPval<- rep(pval.1, c(diff(bp), 1))
    AllPval
    }


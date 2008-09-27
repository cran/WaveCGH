`PlotLatticeOneModel` <-
function(x, Chromosome, PlotChrom=Unique(Chromosome), ModelSelection=FALSE, nsim=200, threshold=0.05, qvalue=FALSE, Panel=c(1, 1), ...)
{
    #
    # Trellis plot for getting break points separtely for each chromosome: 
    #
    #
    # Here we fit one model by considering the residuals obtained from all chromosomes as a single series
    # in our bootstrapping step of obtaining pvalues. 
    # This function will plot the means for different breaks. We can specify whether to smooth 
    # the data before the analysis. If we have only one chromosome number or we want to analyze 
    # an array without specific chromosome, then simply provide "Chromosome= rep(1, length(x))"
    # to get the result. 
    #
    # x           : Array-CGH data
    # Chromosome  : Sequence of chromosomes, eg. 1, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3 and so on.
    # PlotChrom   : Vector of Chromosome numbers which are to be plotted. This is needed because we
    #               fit a single model for all values and then plot individual chromosomes.
    # threshold   : The boundary beyond which we call the region loss or gain.
    # Panel       : How the layout of the lattice plot should we desire.
    # 
        Means <- BreakMeansChr(x = x, Chromosome = Chromosome)
        Pvalues <- PvalBootOneModel(z=x, Chromosome = Chromosome, ModelSelection= ModelSelection, nsim=nsim)
        n<- length(x)
    ind <- 1:n
    br.p <- ind[c(1, diff(Means)) != 0]
      br.p <- Unique(c(br.p, ind[length(x)]))
      bp.n <- as.numeric(na.exclude(br.p))
    Unique.pval<- Pvalues[bp.n[1:(length(bp.n)-1)]]
    if(qvalue){
        Unique.qval<- qvalue(Unique.pval)$qvalues
        Unique.qval.1<- c(Unique.qval, Unique.qval[length(Unique.qval)])
        qValues<- rep(Unique.qval.1, c(diff(bp.n),1))
        xcCode <- rep(7, length(Means))
        xcCode[qValues> threshold & x > 0] <- 5
        xcCode[qValues> threshold & x <= 0] <- 7
        xcCode[qValues< threshold & Means > 0] <- 2
        xcCode[qValues< threshold & Means < 0] <- 3
    }

    if(!qvalue){
        TF<- BenjaminiHochburg(Unique.pval, threshold)
        TF.1<- c(TF, TF[length(TF)])
        TF.th<- rep(TF.1, c(diff(bp.n),1))
        xcCode <- rep(7, length(Means))
        xcCode[TF.th==0 & x > 0] <- 5
        xcCode[TF.th==0  & x <= 0] <- 7
        xcCode[TF.th==1  & Means > 0] <- 2
        xcCode[TF.th==1  & Means < 0] <- 3
    }
    xcCode<-xcCode
    ch.n <- table(Chromosome)
    GeneLabel <- unlist(apply(ch.n, 1, function(x) list(1:x)))
    Ch <- ordered(Chromosome) 
    #
    # We make the vector of 'T' or 'F' corresponding to whether the Chromosome is equal to PlotChrom or not. 
    # 
    LenChrom<-length(PlotChrom)
    v<-matrix(rep(NA,(length(x)*LenChrom)),nrow=LenChrom)
    for (i in 1:LenChrom)v[i,]<-Chromosome==PlotChrom[i]

    NewCol<-rep(NA, length(x))
    for(i in 1:length(x))NewCol[i]<-any(v[,i])
 
    g.df <- data.frame(xr = x[NewCol], xm = Means[NewCol], xcCode = xcCode[NewCol], GeneLabel = GeneLabel[NewCol], Ch = Ch[NewCol])
  
    xyplot(xr ~ GeneLabel | Ch, layout = Panel, strip = function(...)
    strip.default(..., style = 5), 
    panel = function(x, y, subscripts)
    {
        u.m <- Unique(g.df$xm[subscripts])
        xi <- x[c(1, diff(g.df$xm[subscripts])) != 0]
        xi <- Unique(c(xi, x[length(x)]))
        xi <- as.numeric(na.exclude(xi))
        yi <- c(u.m, u.m[length(u.m)])
        ic <- g.df$xcCode[subscripts]
        lsegments(x - 1, 0, x - 1, y, col = ic)
        lsegments(x, 0, x, y, col = ic)
        lsegments(x - 1, y, x, y, col = ic)
        lsegments(x - 1, 0, x, 0, col = ic)
        bp.n <- xi
        g.mean <- yi
        ll <- length(bp.n)
        lsegments(bp.n[ll - 1] - 1, g.mean[ll - 1], bp.n[ll], g.mean[ll - 1], col = 16)
        if(ll > 2) {
            for(i in 1:(ll - 2))
                lsegments(bp.n[i] - 1, g.mean[i], bp.n[i + 1] - 1, g.mean[i], col = 16)
        }
        lsegments(0, 0, 0, g.mean[1])
        if(ll > 2) {
            for(i in 2:(ll - 1)) {
                lsegments(bp.n[i] - 1, 0, bp.n[i] - 1, g.mean[i], col = 16)
                lsegments(bp.n[i] - 1, 0, bp.n[i] - 1, g.mean[i - 1], col = 16)
            }
        }
        lsegments(bp.n[ll], 0, bp.n[ll], g.mean[ll], col = 16)
        lsegments(bp.n[1] - 1, 0, bp.n[ll], 0, col = 16)
    }
    , data = g.df, xlab = "Gene number in chromosome", ylab = "Values", ...)
}


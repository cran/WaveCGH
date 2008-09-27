`PlotAllOneModel` <-
function(x, Chromosome, ModelSelection=FALSE, nsim=200, threshold = 0.05, qvalue=FALSE, ...)
{
    #
    # Here we fit one model by considering the residuals obtained from all chromosomes as a single series
    # in our bootstrapping step of obtaining pvalues.  
    #
    # Plotting all chromosome in one plot:
    #
    # This function will plot the means for different breaks. We can specify whether to smooth 
    # the data before the analysis. If we have only one chromosome number or we want to analyze 
    # an array without specific chromosome, then simply provide "Chromosome= rep(1, length(x))"
    # to get the result. 
    #
    # x           : Array-CGH data
    # Chromosome  : Sequence of chromosomes, eg. 1, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3 and so on.
    # threshold   : The boundary beyond which we call the region loss or gain.
    #  
    Means <- BreakMeansChr(x = x, Chromosome = Chromosome)
    Pvalues <- PvalBootOneModel(z=x, Chromosome = Chromosome,ModelSelection=ModelSelection, nsim=nsim)
     u.m <- Unique(Means)
     n<- length(x)
     ind <- 1:n
        
     br.p <- ind[c(1, diff(Means)) != 0]
     br.p <- Unique(c(br.p, ind[length(x)]))
     bp.n <- as.numeric(na.exclude(br.p))
     Unique.pval<- Pvalues[bp.n[1:(length(bp.n)-1)]]
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
         
    g.mean <- c(u.m, u.m[length(u.m)])
    ic <- xcCode
    ll <- length(bp.n)
    i<-1:length(x)
    plot(c(1,length(x)), c(-max(abs(x)),max(abs(x))), type="n", col = ic, xlab="Clone", ylab="Value", ...)
    segments(i - 1, 0, i - 1, x, col = ic)
    segments(i, 0, i, x, col = ic)
    segments(i - 1, x, i, x, col = ic)
    segments(i - 1, 0, i, 0, col = ic)
    
    # Plotting the values
    lines(c(bp.n[(ll - 1)] - 1, bp.n[ll]), c(g.mean[ll - 1], g.mean[ll - 1]), col = 4)
    if(ll > 2) {
        for(i in 1:(ll - 2))
            lines(c(bp.n[i] - 1, bp.n[i + 1] - 1), c(g.mean[i], g.mean[i]), col = 4)
    }
    lines(c(0, 0), c(0, g.mean[1]), lwd = 3)
    if(ll > 2) {
        for(i in 2:(ll - 1)) {
            lines(c(abs(bp.n[i] - 1), abs(bp.n[i]) - 1), c(0, g.mean[i]), lwd = 1)
            lines(c(abs(bp.n[i]) - 1, abs(bp.n[i]) - 1), c(0, g.mean[i - 1]), lwd = 1)
        }
    }
    lines(c(abs(bp.n[ll]), abs(bp.n[ll])), c(0, g.mean[ll]), lwd = 1)
    abline(v = cumsum(table(Chromosome)), col = 4)
  }


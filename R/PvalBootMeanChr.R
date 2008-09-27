`PvalBootMeanChr` <-
function(x, Chromosome, ModelSelection=FALSE, nsim=200){
 # Finds p-values for testing means each chromosome separately.
 # If we have only one chromosome number or we want to analyze an array without 
 # specific chromosome, then simply provide "Chromosome= rep(1, length(x))"
 # to get the result. 
 #
 # x          : CGH data
 # Chromosome : Chromosome sequence
 # ModelSelection : Whether to select a model or just use AR(1) process
 # nsim           : Number of boostrap simulation
 #
 UniqueChr <- Unique(Chromosome)
 TabChr<-table(Chromosome)
 row.names(TabChr)<-NULL
 dat <- as.list(new.env()) # series as list corresponding to each chromosome
 Pval<-as.list(new.env()) # means as list corresponding to each chromosome
 for (i in 1:length(TabChr)){
     dat[[i]]<-x[Chromosome == UniqueChr[i]]
     Pval[[i]]<-PvalBootMean(dat[[i]],ModelSelection=ModelSelection, nsim=nsim)
     }
 unlist(Pval)
}


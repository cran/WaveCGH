`BreakMeansChr` <-
function(x, Chromosome){
# Finds break points for each chromosome separately and then finds the region means.
# If we have only one chromosome number or we want to analyze an array without 
# specific chromosome, then simply provide "Chromosome= rep(1, length(x))"
# to get the result. 
#
# x          : CGH data
# Chromosome : Chromosome sequence
#
 UniqueChr <- Unique(Chromosome)
 TabChr<-table(Chromosome)
 row.names(TabChr)<-NULL
 dat <- as.list(new.env()) # series as list corresponding to each chromosome
 RegMeans<-as.list(new.env()) # means as list corresponding to each chromosome
 for (i in 1:length(TabChr)){
     dat[[i]]<-x[Chromosome == UniqueChr[i]]
     RegMeans[[i]]<-BreakMeans(dat[[i]])
     }
 unlist(RegMeans)
}


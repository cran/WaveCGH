`BenjaminiHochburg` <-
function(p,q=0.05){
#
# This function finds the significant scores by Benjamini-Hochburg multiple 
# testing method. 
#
# p : p-values
# q : the preselected level
#
 p<-as.numeric(p)
 G<-length(p)
 p.or<-sort(p)
 i.g<-(1:G)/G*q
 TF<-(p.or<=i.g)   # this will give whether p(i)<=iq/G
 i.q<-sum(TF)
 i<- 1:length(p)
 i[order(p)[1:i.q]]<-1
 i[order(p)[(i.q+1):length(p)]]<-0
 i
}


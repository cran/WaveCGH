`Unique` <-
function(x){
#
# This function is different from 'unique'.
#  This finds the unique values in the original order.
#
    a<-duplicated(x)
      i<- 1:length(x)
    i[1]<- 1
      v<-i[as.numeric(a)!=1]
    x[as.numeric(v)]
}


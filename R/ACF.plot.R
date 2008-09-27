`ACF.plot` <-
function(x, cl, ...){
#
# This function will find the acf plot for given x and 
# corresponding chromosome number. 
#
# x  : observations
# cl : chromosome number
#  
    b.m<-BreakMeansChr(x, Chromosome=cl)
    acf(x-b.m, ...)
}


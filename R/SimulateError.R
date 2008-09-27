`SimulateError` <-
function(phi, SigE, n){
    #
    # Simulating Error with AR(1) process:
    #
    # phi  : AR(1) coefficient
    # SigE : Defined below
    # n    : Number of observations to be simulated.
    #
    # We consider the model to be z=mu+e, where e follows AR(1) process.
    # e[i]<- phi*e[i-1]+a[i], where a~NID(0,SigA)
    # Var(e)=SigA^2/(1-phi^2)
    # So fixing the error variance, the innovation variance is given as
    # SigA^2= Var(e)*(1-phi^2) 
    #
    SigA<- sqrt(1-phi^2)*SigE 
    a<-rnorm(n, 0, SigA)
    e<-rep(0, n)
    e[1]<- a[1]*SigE
    for (i in 2:n) e[i]<- phi*e[i-1]+a[i]
        e
}


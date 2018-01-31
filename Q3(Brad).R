install.packages("deSolve")
library(deSolve)

SIR.vector.field <- function(t, vars, parms=c(beta=2,gamma=1)) {
  with(as.list(c(parms, vars)), {
    dS <- -beta*S*I # dS/dt
    dI <- beta*S*I-(1/gamma)*I  # dI/dt
   
    vec.fld <- c(dS=dS, dI=dI)
    return(list(vec.fld)) # ode() requires a list
  })
}

draw.soln <- function(ic=c(S=1,I=0), tmax=1, 
                      times=seq(0,tmax,length.out = 500), 
                      func, parms, ...) {
  soln <- ode(y=ic, times, func, parms)
  lines(times, soln[,"I"], col="blue", lwd=5, ...)
}

tmax <- 10 # end time for numerical integration of the ODE
## draw box for plot:
plot(0,0,xlim=c(0,tmax),ylim=c(0,1),
     type="n",xlab="Time (t)",ylab="Prevalence (I)",las=1)
## initial conditions:
I0 <- 0.010
S0 <- 1 - I0
## draw solutions for several values of parameter beta:
betavals <- c(1.2,1.5,1.8,2,3,4)
for (i in 1:length(betavals)) {
  draw.soln(ic=c(S=S0,I=I0), tmax=tmax,
            func=SIR.vector.field,
            parms=c(beta=betavals[i],gamma=4),
            lty=i
            
          
  )
}
legend("topright",legend=(betavals), lty=1:6)

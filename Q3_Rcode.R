library(deSolve)

SIR.vector.field <- function(t, vars, parms=c(R_0=2,gamma=1)) {
  with(as.list(c(parms, vars)), {
    dS <- -R_0*S*I # dS/dt
    dI <- R_0*S*I-I  # dI/dt
    dR <- I  # dR/dt
    vec.fld <- c(dS=dS, dI=dI, dR=dR)
    return(list(vec.fld)) # ode() requires a list
  })
}

draw.soln <- function(ic=c(S=0.999,I=0.001,R=0), tmax=1, 
                      times=seq(0,tmax,length.out = 500), 
                      func, parms, colour="blue", ... ) {
  soln <- ode(y=ic, times=times, func, parms)
  lines(times, soln[,"I"], col=colour, ... )
}


## Part b
## Plot solutions of the SIR model
tmax <- 10 # end time for numerical integration of the ODE
## draw box for plot:
plot(0,0,xlim=c(0,tmax),ylim=c(0,1),
     type="n",xlab="Time (t)",ylab="Prevalence (I)",las=1)
## initial conditions:
I0 <- 0.001
S0 <- 1 - I0
R0 <-  0
## draw solutions for several values of parameter beta:
Rknot_vals <- c(1.5,2,2.5)
for (i in 1:length(betavals)) {
  draw.soln(ic=c(S=S0,I=I0,R=R0), tmax=tmax,
            func=SIR.vector.field,
            parms=c(R_0=Rknot_vals[i],gamma=1),
            lty=i # use a different line style for each solution
  )
}



# Part c
tmax <- 40 # end time for numerical integration of the ODE
## draw box for plot:
plot(0,0,xlim=c(0,1/gamma_inv*tmax),ylim=c(0,0.5),
     type="n",xlab="Time (t)",ylab="Prevalence (I)",las=1)
## initial conditions:
I0 <- 0.001
S0 <- 1 - I0
R0 <-  0
gamma_inv=4
## draw solutions for several values of parameter beta:
Rknot_vals <- c(1.2,1.5,1.8,2,3,4)
for (i in 1:length(Rknot_vals)) {
  draw.soln(ic=c(S=S0,I=I0,R=R0), tmax=tmax,
            times=1/gamma_inv*seq(0,tmax,length.out = 500),
            func=SIR.vector.field,
            parms=c(R_0=Rknot_vals[i],gamma=1), lwd=2,
            colour=rainbow(length(Rknot_vals)+1)[i] # use a different line colour for each solution
  )
}
legend("topright",legend=Rknot_vals,col=rainbow(length(Rknot_vals)+1),lwd=2)



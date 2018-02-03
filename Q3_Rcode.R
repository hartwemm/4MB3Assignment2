### Part a
library(deSolve)


### Part b
SIR.vector.field <- function(t, vars, parms=c(beta=2,gamma=1)) {
  # Writes vector field of these differential equations
  with(as.list(c(parms, vars)), {
    dS <- -1*beta*S*I # dS/dt
    dI <- beta*S*I-gamma*I  # dI/dt
    dR <- gamma*I  # dR/dt
    vec.fld <- c(dS=dS, dI=dI, dR=dR)
    return(list(vec.fld)) # ode() requires a list
  })
}

draw.soln <- function(ic=c(S=0.999,I=0.001,R=0), tmax=1, 
                      times=seq(0,tmax,length.out = 500), 
                      func, parms=c(R_0=2,gamma_inv=1), colour="blue", ... ) {
  # Solve ode from vector field using deSolve package and 
  # initial condition and parameters given
  with(as.list(parms),{
    gamma=1/gamma_inv
    beta=gamma*R_0
    # times=times/gamma_inv # Convert time to "natural units" of periods of gamma
    soln <- ode(y=ic, times=times, func, parms=c(beta=beta,gamma=gamma)) # Solve the ode 
    lines(times, soln[,"I"], col=colour, ... ) # Add this line to our plot
    return(data.frame(time=times,y=soln[,"I"]))
  })
}


tmax <- 10 # end time for numerical integration of the ODE
## draw box for plot:
plot(0,0,xlim=c(0,tmax),ylim=c(0,1),
     type="n",xlab="Time (t)",ylab="Prevalence (I)",las=1)
## initial conditions:
I0 <- 0.001
S0 <- 1 - I0
R0 <-  0
## Set parameter values
Rknot_vals <- 2.5
gamma_inv <- 1
# Solve ode for these parameters and initial conditions
draw.soln(ic=c(S=S0,I=I0,R=R0), tmax=tmax,
          func=SIR.vector.field,
          parms=c(R_0=Rknot_vals,gamma_inv=gamma_inv,lwd=5)
)


# Part c
tmax <- 40 # end time for numerical integration of the ODE
## initial conditions:
I0 <- 0.001
S0 <- 1 - I0
R0 <-  0
gamma_inv=4
Rknot_vals <- c(1.2,1.5,1.8,2,3,4)
## draw box for plot:
plot(0,0,xlim=c(0,tmax),ylim=c(0,0.5),
     type="n",xlab="Time (t) in days",ylab="Prevalence (I)",las=1)
## draw solution for each value of Rknot:
for (i in 1:length(Rknot_vals)) {
  draw.soln(ic=c(S=S0,I=I0,R=R0), tmax=tmax,
            times=seq(0,tmax,length.out = 500),
            func=SIR.vector.field,
            parms=c(R_0=Rknot_vals[i],gamma_inv=gamma_inv), lwd=2,
            colour=rainbow(length(Rknot_vals)+1)[i] # use a different line colour for each solution
  )
}
legend("topright",legend=Rknot_vals,col=rainbow(length(Rknot_vals)+1),lwd=2)


## Part d
tmax <- 122 # end time for numerical integration of the ODE
euc.dist <- function(x1, x2) sqrt(sum((x1 - x2) ^ 2))
Errors <- matrix(NA,ncol=4,nrow=4*6*6)
y1 <- data.frame(time=philadata$t,y=philadata$pim)
pop=c(900,2000,2500,3000,5000)
for (p in 1:length(pop)) {
  ## initial conditions:
  I0 <- 0.00002
  S0 <- 1 - I0
  R0 <-  0.1
  plot(philadata$t,philadata$pim/5100,pch=20,col="blue",cex=1.2)
  gamma_inv <- 4
  Rknot_vals <- c(1.97,2,2.01,2.03,2.05)
  ## draw solution for each value of Rknot:
  for (i in 1:length(Rknot_vals)) {
    for (j in 1:length(gamma_inv)) {
      y2 <- draw.soln(ic=c(S=S0,I=I0,R=R0), tmax=tmax,
                times=seq(0,tmax,length.out = 488+1),
                func=SIR.vector.field,
                parms=c(R_0=Rknot_vals[i],gamma_inv=gamma_inv[j]), lwd=2, lty=j,
                colour=rainbow(length(Rknot_vals)+1)[i] # use a different line colour for each solution
      )
      func <- y2[which(y2$time %in% y1$time),"y"]
      E <- euc.dist(y1[,2],func)
      Errors[((p-1)*length(Rknot_vals)*length(gamma_inv)+(i-1)*length(gamma_inv)+j),1:4] <- c(pop[p],Rknot_vals[i],gamma_inv[j],E)
    }
  }
  
  legend("topright",legend=Rknot_vals,col=rainbow(length(Rknot_vals)+1),lwd=2)
  title(paste("Population of",pop[p]))
}

plot(NULL)

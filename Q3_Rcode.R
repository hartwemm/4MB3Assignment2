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
tikz("Q3Partc.tex")
tmax <- 40 # end time for numerical integration of the ODE
## initial conditions:
I0 <- 0.001
S0 <- 1 - I0
R0 <-  0
gamma_inv=4
Rknot_vals <- c(1.2,1.5,1.8,2,3,4)
## draw box for plot:
mymar <- par("mar")
mymar["left"] <- mymar["left"] * 0.5
mymar["right"] <- mymar["right"]*0.25
mymar["top"] <- mymar["top"]*0.15
mymar["bottom"] <- mymar["bottom"]*0.5
par(mar=mymar)
plot(0,0,xlim=c(0,tmax),ylim=c(0,0.5),
     type="n",las=1,
     xlab="Time (t)",
     ylab="Prevalence (I)",
     main="SIR curves with varied $\\R_0$ Values")
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
dev.off()

## Part d
tmax <- 122 # end time for numerical integration of the ODE
philadata$t <- seq(1,nrow(philadata))
y1 <- data.frame(time=philadata$t,y=philadata$pim)
euc.dist <- function(x1, x2) sqrt(sum((x1 - x2) ^ 2))
# Errors <- matrix(NA,ncol=4,nrow=4*6*6)


## initial conditions:
I0 <- 0.0001
R0 <-  0.1
S0 <- 1 - I0 - R0
pop <- 8200
plot(philadata$t,philadata$pim/pop,pch=20,col="blue",cex=1.2)
gamma_inv <- 3.5
Rknot_vals <- 1.97

cols <- rainbow(6)
plot(0,0,xlim=c(0,tmax),ylim=c(0,0.1),xlab="Time (t)",
     ylab="Prevalence (I)")
## draw solution for each value of Rknot:
for (i in 1:nrow(Errors)) {
  ### Lowest errors in the order i= 4,1,2,3,5 ####
  I0 <- Errors[i,1][[1]]
  R0 <-  Errors[i,3][[1]]
  S0 <- 1 - I0 - R0
  pop <- Errors[i,4][[1]]
  points(philadata$t,philadata$pim/pop,pch=20,col=cols[i],cex=1.2)
  Rknot_vals <- Errors[i,5]
  gamma_inv <- Errors[i,6]

  
  y2 <- draw.soln(ic=c(S=S0,I=I0,R=R0), tmax=tmax,
                  times=seq(0,tmax,length.out = 488+1),
                  func=SIR.vector.field,
                  parms=c(R_0=Rknot_vals[[1]],gamma_inv=gamma_inv[[1]]), lwd=1,
                  colour=cols[i] # use a different line colour for each solution
  )
  func <- y2[which(y2$time %in% y1$time),"y"]
  E <- euc.dist(y1[,2],func)
  EEE <- c(I0,S0,R0,pop,Rknot_vals,gamma_inv,E)
  EEE
}
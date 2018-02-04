datafile <- "pim_us_phila_city_1918_dy.csv"
philadata <- read.csv(datafile)
philadata$date <- as.Date(philadata$date)

plot(philadata,xlab="Date in 1918",ylab="Prevalence (I)")

lnPI <- log(philadata$pim)
plot(philadata$date,lnPI, pch=21, bg="red", type="p",bty="L",
     ann=FALSE,las=1,xaxs="i",cex=0.65,lwd=0.45,
     xlab="Date in 1918",ylab="ln( Prevalence)")  ## Add red pts

# plot(philadata$date[23:35],lnPI[23:35],xlab="Date in 1918",ylab="ln(P&I)")
m1 <- lm(lnPI[23:35]~philadata$date[23:35])
abline(m1,col="blue",lwd=2)
summary(m1)

########################
## Math 4MB3 Assignment 2, Question 1 - Courtney ##
###############################

#(b) Read in the data
datafile <- "pim_us_phila_city_1918_dy.csv"
philadata <- read.csv(datafile)
philadata$date <- as.Date(philadata$date)

#(c) Plot the data
plot(philadata, pch=21, col="grey", type="l",bty="L",ann=FALSE,las=1,xaxs="i")
points(philadata, pch=21, bg="red",cex=0.65,lwd=0.45)
mtext("Date", side=1, adj=1, line=1.5, font=1, cex=1.5, col='blue') 
mtext("P&I Deaths", side=2, padj=-8, line=-2.5, font=1, las=1, cex=1.5, col='blue')

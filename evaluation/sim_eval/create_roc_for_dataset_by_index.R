#options(echo=TRUE) # if you want see commands in output file
args <- commandArgs(trailingOnly = TRUE)

NUMBER <- args[1]
CMDLINE <- args[2]

data <- read.table("tpr_table.txt", header=FALSE)
#data

#plot(data)

#library(ggplot2)



#max1=max(data$V5)
#max2=max(data$V4)
#max3=max(data$V3)
#max4=max(data$V2)
#maxy=max(max1, max2, max3, max4)


#subset <- t(data.frame(data$V5, data$V4, data$V3, data$V2))
#barplot(subset, legend = c("SweepFinder2", "SweeD", "OmegaPlus", "RAiSD"), col=c("brown3", "chocolate2", "darkgoldenrod3","cadetblue3"), names.arg=data$V1, beside=TRUE, ylab = "Success rate (%)", xlab="Dataset #", log="y")
#title(main = "Success rate for reported location\nin the proximity of the selection target")


# Create Line Chart

# convert factor to numeric for convenience
ntrees <- 4

# get the range for the x and y axis
xrange <- range(data$V1)
yrange <- range(data$V1)

# set up the plot

col=c("brown3", "chocolate2", "darkgoldenrod3","cadetblue3")
col=rev(col)
plot(xrange, yrange, type="l", xlab="FPR",
   ylab="TPR", col="gray", lwd=1)

#colors <- rainbow(ntrees)
##linetype <- c(1:ntrees)

plotchar <- seq(18,18+ntrees,1)
lwd=2
# add lines
#for (i in 1:ntrees) {
  ##tree <- subset(Orange, Tree==1)
lines(data$V1, data$V2, type="l", lwd=lwd,
    lty=1, col=col[1])
lines(data$V1, data$V3, type="l", lwd=lwd,
    lty=1, col=col[2])
lines(data$V1, data$V4, type="l", lwd=lwd,
    lty=1, col=col[3])
lines(data$V1, data$V5, type="l", lwd=lwd,
    lty=1, col=col[4])
#}

# add a title and subtitle
title(paste("ROC curve analysis for Dataset", NUMBER,"\n[", CMDLINE , "]"))

# add a legend
legend("bottomright", legend=c("RAiSD", "OmegaPlus", "SweeD", "SweepFinder2" ), col=col, lty=1)



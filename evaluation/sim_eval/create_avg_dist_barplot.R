args <- commandArgs(trailingOnly = TRUE)

nn <- length(args)

data <- read.table("distance_table.txt", header=FALSE)

max1=max(data$V5)
max2=max(data$V4)
max3=max(data$V3)
max4=max(data$V2)
f = 1.3
maxy=max(max1*f, max2*f, max3*f, max4*f)


subset <- t(data.frame(data$V5, data$V4, data$V3, data$V2))
barplot(subset, legend = c("SweepFinder2", "SweeD", "OmegaPlus", "RAiSD"), col=c("brown3", "chocolate2", "darkgoldenrod3","cadetblue3"), names.arg=args, beside=TRUE, ylab = "Distance (%)", xlab="Dataset #", ylim=c(0,maxy))
title(main = "Average distance between reported location and selection target")


#counts <- table(a)
#barplot(counts, main="Car Distribution by Gears and VS",
#  xlab="Number of Gears", col=c("darkblue","red"), beside=TRUE)




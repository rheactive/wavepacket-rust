library(ggplot2)



for(j in 180:200){

No <- 343*j;

flnm <- paste("results/wavepacket_snapshot-",No,".dat",sep="")

d <- read.table(flnm)

if(j == 0){
#Column height
len <- length(d$V1)
mx <- max(abs(d$V3),abs(d$V4))
}

#Create a new data frame with "long form" like this:

#___x1___y1___"y"
#___x2___y2___"y"
#___x3___y3___"y"
#___x4___y4___"y"
#___x5___y5___"y"
#___x1___z1___"z"
#___x2___z2___"z"
#___x3___z3___"z"
#___x4___z4___"z"
#___x5___z5___"z"

#This is how ggplot2 works

d1 <- data.frame(x = c(d$V1, d$V1, d$V1), f = c(40*d$V2, d$V3, d$V4), group = c(rep("Potential",len),rep("Re(Psi)",len),rep("Im(Psi)",len)))

#Create the plot
g <- ggplot(d1, aes(x = x, y = f, colour = group)) + ylim(-1.0*mx,1.0*mx) +
#Draw lines
geom_line(size=1.5) + 
#Draw points
#geom_point(size = 2.5) +
#Plot style parameters
theme_bw(base_size = 28) + 
theme(aspect.ratio = 0.618, axis.text.x = element_text(colour = "black"), axis.text.y = element_text(colour = "black", angle = 90, hjust = 0.5)) + 
#Legend parameters
theme(legend.title = element_blank(), legend.text = element_text(size = 28), legend.background = element_rect(colour = "black", fill = "white", size = 0.7)) + 
theme(legend.position = "right", legend.spacing.y = unit(10, "pt"), legend.margin=margin(t=0,l=20,b=15,r=20, unit='pt')) + 
guides(colour = guide_legend(byrow = TRUE)) +
#Axes labels
labs(x = "Distance x, nm", y = expression(paste("Wave function")))

#Show the plot on screen
g

#Export the plot

flnm <- paste("plots/wavepacket_snapshot-",No,".png",sep="")

ggsave(flnm, scale = 1, width = 16, height = 10)

}
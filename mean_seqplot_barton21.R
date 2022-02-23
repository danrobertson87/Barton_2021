args <- commandArgs(trailing=TRUE)
library(tools)


library(seqplots)
features<-c(args[1])
tracks<-c(args[2],args[3],args[4],args[5])
a<-getPlotSetArray(refgenome="sacCer3",rm0 = F,ignore_strand = F,add_heatmap = F,
                        tracks=tracks,
                        features=features,
                        bin=50,
                        type="mf",
                        xmin=5000,xmax=5000)

cv = c(args[6], args[7], args[8], args[9])
out = args[10]

pdf(paste0(out, ".pdf")) 
plotAverage(a, plotScale="log2", legend = FALSE, colvec=cv, xlab="Coordinate", ylab="log2(mean calibrated reads)" ,ylim=range(4,12))
dev.off()

png(paste0(out,".png")) 
plotAverage(a, plotScale="log2", colvec=cv, xlab="Coordinate", ylab="log2(mean calibrated reads)" , ylim=range(4,12))
dev.off()

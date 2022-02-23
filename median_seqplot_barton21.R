library(tools)
library(seqplots)
features <- c("sites.bed")
tracks =  c("S1.bw",
			"S2.bw")
			
		                 
options(scipen=9)

a1<-getPlotSetArray(refgenome="sacCer3",rm0 = F,ignore_strand = F,add_heatmap = F,
                                            tracks=tracks, features=features, stat="median",
                                            bin=50,type="mf",xmin=5000,xmax=5000)
median_S1 = a1$data$`sites`$`S1`$mean
median_S2 = a1$data$`sites`$`S2`$mean
x          = a1$data$`sites`$`S2`$all_ind


df = data.frame("Coordinates" = x, "S1" = median_S1, "S2" = median_S2)

library(ggplot2)
library(reshape2)

df.long = melt(df, id.vars="Coordinates")
names(df.long) = c("Coordinate", "Sample", "MedianReadCount")
p = ggplot(df.long, aes(Coordinate, log2(MedianReadCount), colour=Sample)) + geom_line() +
    ylab("log2(median calibrated reads)") + 
    ylim(4, 11) +
    scale_colour_manual(values=c("red", "darkred")) + theme_bw() +
    theme(text = element_text(size=12))

ggsave("median_5kb.pdf", width=15)

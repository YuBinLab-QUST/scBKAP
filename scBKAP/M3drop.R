library(M3Drop)

b <- read.csv('ting_auto.csv',header = F)
Normalized_data <- M3DropCleanData(b, 
                                   is.counts=TRUE, 
                                   min_detected_genes=112)
#the min_detected_genes in each dataset is different,because of different expression values
c <- Normalized_data$data
c <- t(c)
write.csv(c,'ting_m3.csv',row.names = F, col.names = F)

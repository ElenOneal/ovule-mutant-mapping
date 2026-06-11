library(dplyr)

#This script fill create a distribution of heterozygosity in a mapping population
#Excess heterozygosity is probable well contamination
#too little heterozygosity may be a selfed backcross parent
#'AB' is heteryzygous genotype

draw_het_distribution <- function(df){
      hD <- data.frame(colSums(df=='AB'))
      colnames(hD) <- 'nhets'
      hD$fraction <- hD$nhets/nrow(df)
      hD <- hD[-(1:5),]
      hD$sample <- rownames(hD)
      hist(hD$fraction,breaks=100)
      return(hD)
}

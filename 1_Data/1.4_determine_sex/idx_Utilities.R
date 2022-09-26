#March 5 2020
#Written by Else Mikkelsen
#This script parses the results of Samtools idxstats or other .txt files which contain chromosome names in column 1, chromosome length in column 2, and number of mapped sequencing reads in column 3
#Usage:
#Rscript Process_idxstats.R {name_of_idxstats_file.txt} {(optional) length of each sequencing read}

process_idx <- function(filename) {
  idxdata <- read.table(filename)
  colnames(idxdata) <- c("Chromosome", "Length", "Mapped_Reads", "Unmapped_Reads")
  return(idxdata)
}

stats_idx <- function(idxdata, readlength) {
  reportstats <- NULL
  reportstats$Number_Scaffolds <- nrow(idxdata)
  reportstats$Number_Reads <- sum(idxdata$Mapped_Reads)
  reportstats$Genome <- sum(idxdata$Length)
  reportstats$Average_Reads_perscaff <- sum(idxdata$Mapped_Reads/nrow(idxdata))
  reportstats$SD_Reads_perscaff <- sd(idxdata$Mapped_Reads/nrow(idxdata))
  reportstats$seq_depth <- (readlength*reportstats$Number_Reads)/reportstats$Genome
  return(reportstats)
}

plot.depth <- function(idxdata, linear.model) { 
  depthplot <- ggplot(data = idxdata, aes(x = Length, y = Mapped_Reads)) + #create a plot
    geom_point(size = 2, color = "#8a2be2", fill = "#008080", shape = 21) +
    theme_classic() + #remove distracting plot features
    labs(x = "Chromosome Length (bp)", y = "Mapped Sequencing Reads",
         title = "Number of sequencing reads for chromosomes of different lengths") + #add labels that are more meaningful
    theme(plot.title = element_text(hjust = 0.5))+  #make title centred
    geom_abline(slope = linear.model$coefficients[2], intercept = linear.model$coefficients[1])
  return(depthplot)
}

#This function calculates the correlation between wine density and residual sugar content
depth.cor <- function(idxdata) {
  correlation <- cor(idxdata$Length, idxdata$Mapped_Reads)
  return(correlation)
}

#This function calculates the covariance between wine density and residual sugar content
depth.cov <- function(idxdata) {
  covariance <- var(idxdata$Length, idxdata$Mapped_Reads)
  return(covariance)
}

#This function creates a linear model between wine density and residual sugar content
depthlm <- function(idxdata) {
  lin.model <- lm(formula = Mapped_Reads ~ Length, data = idxdata) #create linear model
  return(lin.model) #return linear model
}

guessZ <- function(idxdata, linmod) {
  points <- linmod$coefficients[2]*idxdata$Length + linmod$coefficients[1]
  #In females, the Z (or W) chromosome should show the greatest deviation from the expected value of sequencing depth, and it should always be lower than expected
  deviation <- points - idxdata$Mapped_Reads
  names(deviation) <- idxdata$Chromosome
  putativeZ <- names(deviation)[deviation == max(deviation)]
  return(putativeZ)
}


outlierscaffs <- function(idxdata, linmod) {
  #Make sure that this is a statistically significant outlier
  resids <- residuals(linmod)
  names(resids) <-idxdata$Chromosome
  outliers <- names(abs(resids))[abs(resids)>4*sd(resids)&abs(as.numeric(resids))>0.2*idxdata$Mapped_Reads]
  return(outliers)
}

Zbinomial <- function(idxdata, putativeZ, linmod) {
  #Select data for the Z chromosome
  Zdata <- idxdata[idxdata$Chromosome==putativeZ,]
  #Create expected data based on lengths of chromosomes
  points <- linmod$coefficients[2]*idxdata$Length + linmod$coefficients[1]
  #This value is expected to be equal to 0.5:
  #(Zdata$Mapped_Reads/Zdata$Length)/(points[idxdata$Chromosome==putativeZ]/Zdata$Length)
  #We will test whether it deviates significantly
  #To assess significance, we should test if there is a significant deviation from 0.5
  #This is a test of whether, given the sample size, the (# of mapped reads)/(# of expected reads) is significantly different from 0.5
  #This is essentially a binomial test. If we consider a mapped read to be a "success" and the number of expected reads for a chromosome of its size to be the "number of trials" and the expected proportion of "successes" to be 50%:
  binomial <- binom.test(x=Zdata$Mapped_Reads, n=as.integer(points[idxdata$Chromosome==putativeZ]), p=0.5)
  return(binomial)
}


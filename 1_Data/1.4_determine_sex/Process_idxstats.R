#March 5 2020
#Written by Else Mikkelsen
#This script parses the results of Samtools idxstats or other .txt files which contain chromosome names in column 1, chromosome length in column 2, and number of mapped sequencing reads in column 3
#Usage:
#Rscript Process_idxstats.R {name_of_idxstats_file.txt} {(optional) length of each sequencing read}

#Load environment
library(ggplot2)
source("idx_Utilities.R")

#Parse command line arguments
command <- commandArgs(trailingOnly = T)

#end program if there are no command line arguments given and print an error message
if (length(command)==0){
  cat("Error: This scripts requires the path to the idxstats file as an argument\n")
  cat("Usage:\n")
  cat("Rscript Process_idxstats.R {name_of_idxstats_file.txt} {(optional) length of each sequencing read}\n")
  quit()
}

#end program if there are too many command line arguments given and print an error message.
#If the optional second argument is given (length of sequencing reads, measured in basepairs), 
#then assign this value to a variable
#Otherwise, assume a common default value
if (length(command)>1){ #check that there are multiple command line arguments
  command2 <- as.numeric(command[2]) #make the second argument numeric
  readlength <- command2 #assign read length to a variable
  if (!is.numeric(command2)|is.na(command2)){ #check that it was successfully converted to numeric
    cat("Error: second command line argument is not an integer. This script can only process one file at a time.\n")
    cat("The second command line argument should be the integer length of the sequencing reads (optional)\n")
    quit()
  }
} else {
  readlength <- 150 #if no second argument given, assume a default read length
}

#report the read length that will be used
cat("Using a read length of", readlength, "bp.\nAdjust second command line argument if this is not correct.\n")

#end program if there are too many command line arguments given and print an error message
if (length(command)>2){
  cat("Error: This script can only process one file at a time.\n")
  cat("The second command line argument should be the integer length of the sequencing reads (optional)\n")
  quit()
}


#####################
idxstats <- command[1] #get the filename to be used

cat("The file being processed is", idxstats, "\n") #report the sample being analyzed

#First, load the data
idxdata <- process_idx(idxstats)

cat("\n-----------------------------\n")

#check integrity of data - should not be any unmapped reads in the file, but if there are, the program should not quit
if (sum(idxdata$Unmapped_Reads)>0){
  cat("Warning: Did you know that there are unmapped reads in this data?\n")
}
if (nrow(idxdata)==0|!exists("idxdata")){
  cat("Error: This file appears to be empty\n")
  quit() #exit the program if the file could not be read
}
if (sum(idxdata$Mapped_Reads)==0){
  cat("Error: There appear to be no mapped reads in this file\n")
  quit()
}

#report basic stats about the data: number of scaffolds, length of genome, number of reads, sequencing depth, and the average and standard deviation of number of reads per chromosome
reportstats <- stats_idx(idxdata, readlength)
cat("There are", reportstats$Number_Scaffolds, "scaffolds in this dataset, with a total length of", reportstats$Genome, "bp\n")
cat("There are", reportstats$Number_Reads, "mapped reads for this individual, for an average sequencing depth of", reportstats$seq_depth, "X \n")
cat("The mean number of reads per chromosome is", reportstats$Average_Reads, "with a standard deviation of", reportstats$SD_Reads_perscaff, "\n")

#report the correlation and covariance between scaffold length and mapped reads, without removing the Z chromosome
cat("The correlation between scaffold length and mapped reads is", depth.cor(idxdata), "\n")
cat("The covariance between scaffold length and mapped reads is", depth.cov(idxdata), "\n")
cat("\n-----------------------------\n\n")


#create a linear model of the relationship between scaffold length and mapped reads.
linmod <- depthlm(idxdata)
cat("Created a linear model of the relationship between scaffold length and mapped reads:\n")
print(summary(linmod))

#Visualize the data, and plot a line based on the linear model. 
length_vs_reads <- plot.depth(idxdata, linmod)
pdf(paste(idxstats, "plots.pdf", sep = "")) #Save the plot to a pdf.
print(length_vs_reads)

#Guess which chromosome is the Z chromosome, and detect which scaffolds are outliers based on the linear model
putativeZ <- guessZ(idxdata, linmod)
outliers <- outlierscaffs(idxdata, linmod)

#Report the Z chromosome, if the putative Z chromosome is in the set of outliers, then report it to the user. Otherwise, tell them it is probably a male bird
if (putativeZ %in% outliers){
  foundZ <- TRUE
  cat("Congratulations! It looks like the Z chromosome is", putativeZ, "\n")
  if (length(outliers)>1){ #report other outlier scaffolds so that the user can investigate them if needed
    cat("You might also want to investigate all these scaffolds:", outliers[outliers != putativeZ], "\n")
  }
} else {
  foundZ <- FALSE
  cat("I could not find the Z chromosome. My closest guess is", putativeZ, "\n")
  cat("This is probably a male bird (or female mammal)\n")
  if (length(outliers)>0){ #Even though we did not detect the Z chromosome, report any outliers to the user so they can investigate them if needed
    cat("You might want to investigate these scaffolds:", outliers, " because they appear to be outliers\n")
  }
}

cat("\n-----------------------------\n")

#If we found the Z chromosome, report the stats without the Z chromosome
if (foundZ){
  autosomedata <- idxdata[idxdata$Chromosome !=putativeZ,] #create a dataset with the Z removed
  #report the correlation and covariance between scaffold length and mapped reads, after removing the Z chromosome
  cat("The correlation between scaffold length and mapped reads is", depth.cor(autosomedata), "after removing the Z\n")
  cat("The covariance between scaffold length and mapped reads is", depth.cov(autosomedata), "after removing the Z\n")
  #create a linear model of the relationship between scaffold length and mapped reads, without the Z.
  autolinmod <- depthlm(autosomedata)

  #Visualize the data, and plot a line based on the linear model
  cat("A plot of sequencing reads vs chromosome length has been produced (with and without the putative Z, if detected.)\n")
  cat("You may want to inspect it visually to confirm that it follows the expected pattern\n")
  plot.no.Z <- plot.depth(autosomedata, autolinmod)
  print(plot.no.Z) #save plot to the pdf file
  
  #test whether the relative sequencing depth deviates from the expected value
  binomial <- Zbinomial(idxdata, putativeZ, autolinmod)  #do a binomial test
  if (binomial$p.value<0.05){ #check if the result was significant
    cat("\n-----------------------------\n")
    cat("Warning: the depth of coverage deviates significantly from the expected value of half the autosomal sequencing depth\n")
    cat("Binomial test: p=", binomial$p.value, "\n") #report the p value
    if (binomial$estimate>0.5){ #check whether the estimate was higher or lower than expected
      cat("The sequencing depth of the Z is", binomial$estimate, "times the depth of the other chromosomes\n")
      cat("Warning: a very elevated value here could indicate inaccurate read mapping\n")
    }
    if (binomial$estimate<0.5){
      cat("The sequencing depth of the Z is", binomial$estimate, "the depth of the other chromosomes\n")
    }
  }
  if (binomial$p.value>0.05){ #report to the user if the result of the binomial test was not significant
    cat("The relative sequencing depth of the Z (", binomial$estimate, ") does not deviate significantly from expectation (0.5)\n" )
    cat("Binomial test: p=", binomial$p.value, "\n")
  }
}  

#close pdf
dev.off()

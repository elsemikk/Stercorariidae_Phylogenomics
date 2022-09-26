#This file is to be used after Process_idxstats.R to confirm the sex of male birds after the Z chromosome has been identified
#Else Mikkelsen Mar 6 2020
#Load environment
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

idxstats <- command[1] #get the filename to be used

cat("The file being processed is", idxstats, "\n") #report the sample being analyzed

#First, load the data
idxdata <- process_idx(idxstats)

cat("\n-----------------------------\n")

putativeZ <- command[2]

autosomedata <- idxdata[idxdata$Chromosome !=putativeZ,] #create a dataset with the Z removed
#create a linear model of the relationship between scaffold length and mapped reads, without the Z.
autolinmod <- depthlm(autosomedata)

#test whether the relative sequencing depth deviates from the expected value
binomial <- Zbinomial(idxdata, putativeZ, autolinmod)  #do a binomial test
if (binomial$p.value<0.05){ #check if the result was significant
  cat("\n-----------------------------\n")
  cat("Binomial test is significant: p=", binomial$p.value, "\n") #report the p value
  if (binomial$estimate>0.7){ #check whether the estimate was higher or lower than expected
    cat("The sequencing depth of the Z is", binomial$estimate, "times the depth of the other chromosomes\n")
    cat("Yes, this appears to be a male\n")
  }
  if (binomial$estimate<0.7){
    cat("The sequencing depth of the Z is", binomial$estimate, "the depth of the other chromosomes\n")
    cat("I think this is actually a female\n")
    }
}
if (binomial$p.value>0.05){ #report to the user if the result of the binomial test was not significant
  cat("The relative sequencing depth of the Z (", binomial$estimate, ") does not deviate significantly from expectation (0.5)\n" )
  cat("Binomial test: p=", binomial$p.value, "\n")
  cat("This is probably a female")
}

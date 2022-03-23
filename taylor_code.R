
# INPUT IS TEXT FILE
get_phage_bact_coverages <- function(name_of_file) {
  testPhage <- read.delim(name_of_file, header=FALSE, sep=":")
  bactMatches <- strsplit(testPhage[1,2],", ")
  bactCovs <- as.numeric(unlist(strsplit(testPhage[2,2],",")))
  phagePreds <- strsplit(testPhage[3,2],", ")
  phageCovs <- as.numeric(unlist(strsplit(testPhage[4,2],",")))
  my_list <- list("bact"=bactMatches,"bcovs"=bactCovs,"phage"=phagePreds,"pcovs"=phageCovs)
  return(my_list)
}

bootstrapping <- function(name_of_file,name_for_plot,threshold_value){
  b_p_sample <- get_phage_bact_coverages(name_of_file)
  bootsampleBact <-sample(b_p_sample$bcovs, size = 10000, replace = TRUE)
  summary(bootsampleBact)
  d<-density(bootsampleBact)
  plot(d, main = name_for_plot)
  
  threshold<-as.numeric(quantile(bootsampleBact, probs = threshold_value))
  #rebels are outside the threshold value percentile (99th or 95th etc.)
  rebels<-ifelse(b_p_sample$pcovs>threshold, b_p_sample$pcovs, NA)
  rebel_coverages<-rebels[!is.na(rebels)]
  NonNAindex<-which(!is.na(rebels))
  rebel_names<-b_p_sample$phage[[1]][NonNAindex]
  combo_rebels<-cbind(rebel_names,rebel_coverages)
  
  return(combo_rebels)
}
myArgs <- commandArgs(trailingOnly = TRUE)
file_name<-myArgs[1]
plot_name<-myArgs[2]
threshold_number<-as.numeric(myArgs[3])

rebels<-bootstrapping(file_name,plot_name,threshold_number)

write.csv(rebels,"phages_exceeding_threshold.csv", row.names = FALSE)

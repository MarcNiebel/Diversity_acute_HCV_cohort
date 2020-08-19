#Created by Marc Niebel July 2020
#This script combines all the patients genome region
#diversity measurements into respective csv files

#Pulling out the patient number and genotype
library(stringi)
#Combining clinical outcome with diversity measurements
library(dplyr)

#All datafiles listed
data_shannon_files <-list.files(path ="Output",pattern="*shannon_by_region.csv",full.names = TRUE)
data_gini_simpson_files <- list.files(path="Output",pattern="*gini_simpson_by_region.csv",full.names = TRUE)
data_uncorrected_p_dist_files <- list.files(path="Output",pattern="*p_dist_by_region.csv",full.names = TRUE)

#Reading in the clinical outcome for the patients analysed
clinical_outcome <- read.csv("Data/clinical_outcome.csv")

#extracting the patient number
patient_number <- lapply(data_shannon_files,function(name) as.character(stri_extract(name, regex="\\w\\d{1,}")))
#Extracting the genotype
infected_genotype <- lapply(data_shannon_files,function(genotype) as.character(stri_extract(genotype,regex="\\d{1}[a-z]")))
#Reading in all the csv files into a list
data_shannon <- lapply(data_shannon_files,read.csv)
data_gini_simpson <-lapply(data_gini_simpson_files,read.csv)
data_p_dist <- lapply(data_uncorrected_p_dist_files,read.csv)
#Adding each patient number to respective list
add_patient_number_shannon <- mapply(`[<-`,data_shannon,'patient_no',value=patient_number,SIMPLIFY = FALSE)
add_patient_number_gini_simpson <- mapply(`[<-`,data_gini_simpson,'patient_no',value=patient_number,SIMPLIFY = FALSE)
add_patient_number_p_dist <- mapply(`[<-`,data_p_dist,'patient_no',value=patient_number,SIMPLIFY = FALSE)
#Adding each genotype to respective list
add_genotype_shannon <- mapply(`[<-`,add_patient_number_shannon,'infected_gt',value=infected_genotype,SIMPLIFY = FALSE)
add_genotype_gini_simpson <- mapply(`[<-`,add_patient_number_gini_simpson,'infected_gt',value=infected_genotype,SIMPLIFY = FALSE)
add_genotype_p_dist <- mapply(`[<-`,add_patient_number_p_dist,'infected_gt',value=infected_genotype,SIMPLIFY = FALSE)
#Combining the list into one dataframe
shannon_dataframe <-do.call(rbind,add_genotype_shannon)
gini_simpson_dataframe <-do.call(rbind,add_genotype_gini_simpson)
p_dist_dataframe <- do.call(rbind,add_genotype_p_dist)
#Adding clinical outcome
shannon_dataframe_co <-full_join(shannon_dataframe,clinical_outcome)
gini_simpson_dataframe_co <-full_join(gini_simpson_dataframe,clinical_outcome)
p_dist_dataframe_co <- full_join(p_dist_dataframe,clinical_outcome)

#Writing to file
write.csv(shannon_dataframe_co,"Data/shannon_entropy_genome_regions.csv")
write.csv(gini_simpson_dataframe_co,"Data/gini_simpson_genome_regions.csv")
write.csv(p_dist_dataframe_co,"Data/uncorrected_p_dist_genome_regions.csv")

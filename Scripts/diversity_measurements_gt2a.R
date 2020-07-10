#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = TRUE)

#Created by Marc Niebel in June 2020
#Calculates shannon entropy,gini-simpson and uncorrected p-distance for
#alignment windows generated across the HCV genome

#Needed for most of the alignment manipulations
library(QSutils)
#Used to sort files
library(gtools)
#Used to pull out amino acid cooridnates
library(stringi)
#Needed to generate graphs and manipulations
library(ggpubr)
#join dataframes and summarising
library(dplyr)


#The nucleotide alignment locations for a patient 
alignment_files <- list.files(args[1],pattern = "*.goodfna",full.names = TRUE)
#Set them in order
alignment_files<- mixedsort(alignment_files)
#Extracting the amino acid coordinates from each file and putting them in a list
amino_acid_cooridnates <- lapply(alignment_files, function(name) as.character(stri_match_last(name, regex="\\d{1,}\\w_\\d{1,}")))
#Regular expression is not pulling out
#respective windows so being modified
amino_acid_cooridnates[[1]]<-"1_9"
amino_acid_cooridnates[[264]] <-("2400_2400h")
#Making ranges of desired windows
min1 <- as.list(seq(1,469,by=9))
min2 <- as.list("477")
min3 <-as.list(seq(485,566,by=9))
min4 <- as.list("573b")
min5 <- as.list(seq(582,2229,by=9))
min6 <- as.list(seq(2242,2377,by=9))
min7 <- as.list(seq(2382,2391,by=9))
min8 <- as.list("2400")
min9 <-as.list("2400i")
min10 <-as.list("2401a")
min11 <- as.list("2407")
min12 <-as.list(seq(2418,2994,by=9))
max1 <- as.list(seq(9,468,by=9))
max2 <- as.list("476")
max3 <- as.list(seq(484,565,by=9))
max4 <-as.list("573")
max5 <- as.list(seq(581,2228,by=9))
max6 <- as.list(seq(2241,2376,by=9))
max7 <-as.list(seq(2381,2399,by=9))
max8 <- as.list("2400h")
max9 <- as.list("2401")
max10 <- as.list("2406")
max11 <-as.list("2417")
max12 <- as.list(seq(2426,3002,by=9))

#Combining the lists
mins <- c(min1,min2,min3,min4,min5,min6,min7,min8,min9,min10,min11,min12)
maxs <- c(max1,max2,max3,max4,max5,max6,max7,max8,max9,max10,max11,max12)

desired_amino_acid_windows <- as.data.frame(paste0(mins,sep="_",maxs))
names(desired_amino_acid_windows)[1] <- "amino_acid_windows"
#Regions(aa); Core:1:191(1:21),E1:192:383 (22:43),
#E2: 384:746(44:83),P7:747:809(84:90,NS2:810:1026 (91:114)
#NS3:1027:1657(115:185),NS4A:1658:1711(186:191)
#NS4B:1712:1972(192:220),NS5A:1973:2420(221:271),NS5B:2421:3011(272:336)
#Adding a column of numbers to dataframe
desired_amino_acid_windows$numbers <- 1:nrow(desired_amino_acid_windows)
desired_amino_acid_windows <- desired_amino_acid_windows %>% 
    mutate(Genome_regions=case_when(numbers < 22 ~ "Core",
                                    numbers >= 22 & numbers <= 43 ~ "E1",
                                    numbers >=44 & numbers <= 83 ~ "E2",
                                    numbers >=84 & numbers <=90 ~ "P7",
                                    numbers >= 91 & numbers <=114 ~"NS2",
                                    numbers >= 115 & numbers <= 185 ~"NS3",
                                    numbers >= 186 & numbers <= 191 ~"NS4A",
                                    numbers >=192 & numbers <= 220 ~"NS4B",
                                    numbers >=221 & numbers <= 271 ~"NS5A",
                                    TRUE ~ "NS5B"))
#Removal of number column
desired_amino_acid_windows <- desired_amino_acid_windows[,-2]

#Will be needed as a vector later on
amino_acid_cooridnates <- as.character(amino_acid_cooridnates)
#Making a list of sequence alignment files
sequence_files <- lapply(alignment_files,readDNAStringSet)
#Generating the haplotypes for each alignment
reads_collapsed <- lapply(sequence_files,Collapse)
#Pulling out the frequencies
haplotype_freq <- sapply(reads_collapsed,"[",1)
#Calculating the shannon entropy for each window 
shannon_entropy <- sapply(haplotype_freq,Shannon)
#Creating this as a dataframe
shannon_entropy_df <- as.data.frame(shannon_entropy)
#Adding the amino acid windows to the shannon entropy dataframe
shannon_entropy_df$amino_acid_windows <- amino_acid_cooridnates
#Results of desired windows and actual windows
merge_desired_actual_windows_shannon <- left_join(desired_amino_acid_windows,shannon_entropy_df,by="amino_acid_windows")
#All NA are replaced with 0
merge_desired_actual_windows_shannon[is.na(merge_desired_actual_windows_shannon)] <- 0
plot1 <- ggscatter(merge_desired_actual_windows_shannon,x="amino_acid_windows",y="shannon_entropy",
                   xlab="Nucleotide windows across HCV genome", ylab = "Shannon entropy",
                   title="Shannon entropy across HCV open reading frame",color = "#00AFBB")+
                   expand_limits(x=c(-2,336))+rremove("x.text")+rremove("x.ticks")
#Average shannon entropy
average_shannon_entropy <- mean(shannon_entropy)
average_shannon_entropy <- signif(average_shannon_entropy,digits = 2)
average_shannon_entropy <- paste("Average Shannon entropy",average_shannon_entropy)
#Average by genome region
average_shannon_by_region <- merge_desired_actual_windows_shannon %>%
    group_by(Genome_regions) %>%
    dplyr::summarise(mean_by_region=mean(shannon_entropy))

#Gini-Simpson Index(In verse of Simpson Index() which is that 0 is highly diverse.
#Therefore in this instance it is reported as the inverse where 1 is highly diverse)
gini_simpson <-sapply(haplotype_freq,GiniSimpson)
gini_simpson_df <- as.data.frame(gini_simpson)
gini_simpson_df$amino_acid_windows <- amino_acid_cooridnates
#Results of desired windows and actual windows
merge_desired_actual_windows_gini <- left_join(desired_amino_acid_windows,gini_simpson_df,by="amino_acid_windows")
#All NA are replaced with 0
merge_desired_actual_windows_gini[is.na(merge_desired_actual_windows_gini)] <- 0
plot2 <- ggscatter(merge_desired_actual_windows_gini,x="amino_acid_windows",y="gini_simpson",
          xlab="Nucleotide windows across HCV genome", ylab = "Gini Simpson Diversity Index",
          title="Gini Simpson Diveristy Index across HCV open reading frame",color = "#E7B800")+
          expand_limits(x=c(-2,336))+rremove("x.text")+rremove("x.ticks")
#Average Gini-Simpson Index
average_gini_simpson <- mean(gini_simpson)
average_gini_simpson <- signif(average_gini_simpson,digits = 2)
average_gini_simpson <- paste("Average Gini-Simpson Index:",average_gini_simpson) 
#Average by genome region
average_gini_simpson_by_region <- merge_desired_actual_windows_gini %>%
    group_by(Genome_regions) %>%
    dplyr::summarise(mean_by_region=mean(gini_simpson))

#Functional diveristy
alignments_only <- sapply(reads_collapsed,"[",2)
#Raw distances between sequences
dist <- lapply(alignments_only,DNA.dist)
#nucleotide diversity by entity(average difference between haplotypes in the alignment)
nucleotide_diveristy_entity <- sapply(dist, NucleotideDiversity)
nucleotide_diveristy_entity_df <- as.data.frame(nucleotide_diveristy_entity)
nucleotide_diveristy_entity_df$amino_acid_windows <- amino_acid_cooridnates
#Results of desired windows and actual windows
merge_desired_actual_windows_p_dist <- left_join(desired_amino_acid_windows,nucleotide_diveristy_entity_df,by="amino_acid_windows")
#All NA are replaced with 0
merge_desired_actual_windows_p_dist[is.na(merge_desired_actual_windows_p_dist)] <- 0
plot3 <- ggscatter(merge_desired_actual_windows_p_dist,x="amino_acid_windows",y="nucleotide_diveristy_entity",
                   xlab="Nucleotide windows across HCV genome",ylab ="Uncorrected p-distance",
                   title="Uncorrected p-distance across HCV open reading frame",color = "#FC4E07")+
                   expand_limits(x=c(-2,336))+rremove("x.text")+rremove("x.ticks")
#Average uncorrected p-distance
average_uncorrected_p_distance <- mean(nucleotide_diveristy_entity)
average_uncorrected_p_distance <- signif(average_uncorrected_p_distance,digits=2)
average_uncorrected_p_distance <- paste("Average uncorrected p-distance:",average_uncorrected_p_distance)

#Average by genome region
average_uncorrected_p_dist_by_region <- merge_desired_actual_windows_p_dist %>%
    group_by(Genome_regions) %>%
    dplyr::summarise(mean_by_region=mean(nucleotide_diveristy_entity))

#Extract patient
patient_number <- stri_extract(alignment_files[1], regex="\\w\\d{1,}")
#Extract genotype
infected_genotype <- stri_extract(alignment_files[1],regex="\\d{1}[a-z]")
#PDF output(All three graphs on one page)
ggarrange(plot1,plot2,plot3,nrow=3,ncol=1)%>%
  ggexport(filename =paste("Output/",patient_number,"Diversity plots.pdf"))

#Average outputs
sink(paste("Output/",patient_number,"Averages.txt"))
print(average_shannon_entropy)
print(average_gini_simpson)
print(average_uncorrected_p_distance)
sink(file=NULL)

#Average outputs by genome region
write.csv(average_shannon_by_region,paste0("Output/",patient_number,"_",infected_genotype,"_","average_shannon_by_region.csv"))
write.csv(average_gini_simpson_by_region,paste0("Output/",patient_number,"_",infected_genotype,"_","average_gini_simpson_by_region.csv"))
write.csv(average_uncorrected_p_dist_by_region,paste0("Output/",patient_number,"_",infected_genotype,"_","average_uncorrected_p_dist_by_region.csv"))


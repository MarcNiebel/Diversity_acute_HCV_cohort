#Created by Marc Niebel July 2020
#This script visualises by errorplot and carries out statistical
#analysis by wilcoxon test(non-parametric) whether the median is
#different between the spontaneous clears and progressors

#Diversity measurment data
data_shannon <-read.csv("Data/shannon_entropy_genome_regions.csv")
data_gini_simpson <- read.csv("Data/gini_simpson_genome_regions.csv")
data_p_dist <- read.csv("Data/uncorrected_p_dist_genome_regions.csv")
#Change order of genome regions
levels(data_shannon$Genome_regions) <- c("Core","E1","E2","P7","NS2",
                                        "NS3","NS4A","NS4B","NS5A","NS5B")
levels(data_gini_simpson$Genome_regions)<- c("Core","E1","E2","P7","NS2",
                                             "NS3","NS4A","NS4B","NS5A","NS5B")
levels(data_p_dist$Genome_regions)<- c("Core","E1","E2","P7","NS2",
                                       "NS3","NS4A","NS4B","NS5A","NS5B")

#Errorplot with median and interquartile range
shannon_errorplot <- ggerrorplot(data_shannon,x="Genome_regions",y="mean_by_region",
                                 xlab = "Genome Regions",ylab = "Average Shannon Entropy",
                                 legend.title = "Clinical Outcome",
                                 desc_stat = "median_iqr",color="Clinical.Outcome")+
                                 stat_compare_means(aes(group=Clinical.Outcome),label="p.signif")

gini_simpson_errorplot <- ggerrorplot(data_gini_simpson,x="Genome_regions",y="mean_by_region",
                                      xlab = "Genome Regions",ylab = "Average Gini-Simpson Index",
                                      legend.title = "Clinical Outcome",
                                      desc_stat = "median_iqr",color="Clinical.Outcome")+
                                      stat_compare_means(aes(group=Clinical.Outcome),label="p.signif")

p_dist_errorplot <- ggerrorplot(data_p_dist,x="Genome_regions",y="mean_by_region",
                                      xlab = "Genome Regions",ylab = "Average uncorrected p-distance",
                                      legend.title = "Clinical Outcome",
                                      desc_stat = "median_iqr",color="Clinical.Outcome")+
                                      stat_compare_means(aes(group=Clinical.Outcome),label="p.signif")
multi_page_errorplot <- ggarrange(shannon_errorplot,gini_simpson_errorplot,
                                    p_dist_errorplot,nrow=1,ncol=1)
#Output of the three diversity measurements as a whole
ggexport(multi_page_errorplot,filename ="Output/Errorplots_for_diversity_by_genome_regions.pdf")

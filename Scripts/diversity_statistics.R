#Created by Marc Niebel July 2020
#This script allows for visualisation of diversity measurements
#across the genome and associated statistics

#Library required to produce the box plot and statistics
library(ggpubr)
library(dplyr)

#Diversity measurements data
data <- read.csv("./Data/diversity_measurements.csv")

Shannon_entropy_data <- ggboxplot(data,y="Average.Shannon.Entropy",
                                  x="Clinical.Outcome",color="Clinical.Outcome",
                                  xlab = "Clinical Outcome",ylab="Average Shannon Entropy",legend="none",add = "jitter")+
                                  stat_compare_means(label="p.format",hjust=-1.0)
Gini_simpson_data <- ggboxplot(data, y="Average.Gini.Simpson.Index",x="Clinical.Outcome",
                                       color="Clinical.Outcome",
                                       xlab = "Clinical Outcome",ylab="Average Gini-Simpson Index",legend="none",add = "jitter")+
                                       stat_compare_means(label="p.format",hjust=-1.0)
Uncorrected_p_distance <- ggboxplot(data, y="Average.uncorrected.p.distance",x="Clinical.Outcome",
                                       color="Clinical.Outcome",
                                       xlab = "Clinical Outcome",ylab="Average Uncorrected p-distance",legend="none",add = "jitter")+
                                       stat_compare_means(label="p.format",hjust=-1.0)
multi_page <- ggarrange(Shannon_entropy_data,Gini_simpson_data,Uncorrected_p_distance,nrow=1,ncol=1)
#Output of the three diversity measurements as a whole
ggexport(multi_page,filename ="Output/Summary_statistics_for_diversity.pdf")
#Only for genotypes with > 1 result for each genotype ie gt1a and gt3a
gt1a_gt3a <- data %>% filter(Infected.genotype =="1a"|Infected.genotype=="3a")
#Diversity measurements for all three based on genotype
statistics_by_genotype_shannon <-ggboxplot(gt1a_gt3a,y="Average.Shannon.Entropy",x="Clinical.Outcome",
                                                color = "Clinical.Outcome",add = "jitter",
                                                facet.by = "Infected.genotype",legend.title="Clinical Outcome",
                                                ylab="Average Shannon Entropy ")+rremove("x.text")+rremove("x.title")+
                                                stat_compare_means(label = "p.format",hjust=-0.75)
statistics_by_genotype_gini_simpson <- ggboxplot(gt1a_gt3a,y="Average.Gini.Simpson.Index",x="Clinical.Outcome",
                                   color = "Clinical.Outcome",add = "jitter",
                                   facet.by = "Infected.genotype",legend.title="Clinical Outcome",
                                   ylab="Average Gini-Simpson Index")+rremove("x.text")+rremove("x.title")+
                                   stat_compare_means(label = "p.format",hjust=-0.75)
statistics_by_genotype_uncorrected_p_distance <- ggboxplot(gt1a_gt3a,y="Average.uncorrected.p.distance",x="Clinical.Outcome",
                                                 color = "Clinical.Outcome",add = "jitter",
                                                 facet.by = "Infected.genotype",legend.title="Clinical Outcome",
                                                 ylab="Average uncorrected p-distance")+rremove("x.text")+rremove("x.title")+
                                                 stat_compare_means(label = "p.format",hjust=-0.75)
multi_page_by_genotype <- ggarrange(statistics_by_genotype_shannon,statistics_by_genotype_gini_simpson,
                                    statistics_by_genotype_uncorrected_p_distance,nrow=1,ncol=1)
#Output of the three diversity measurements as a whole
ggexport(multi_page_by_genotype,filename ="Output/Summary_statistics_for_diversity_by_genotype.pdf")

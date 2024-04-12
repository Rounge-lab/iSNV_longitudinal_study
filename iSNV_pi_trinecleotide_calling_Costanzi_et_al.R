rm(list = ls())

library("ggplot2")
library("ggpubr")
library("purrr")
library("dplyr")
library("tidyr")
library("tidyverse")
library("readxl")

#define variables
HPV_subtype <-"HPV16"
setwd("/longitudinal_study_HPV16_2023")
cutoff_meancov<-300

#define functions
# Function for finding column name for specific value
maxN <- function(x, reachback = 0){
  # reachback = 0 is maximum, 1 is second to last, 2 is third to last and so on
  len <- length(x)
  if(reachback > len){
    error('You can not overreach the number of variables.')
  }
  names(sort(x, decreasing = TRUE)[1 + reachback])
}

# Function for finding the second largest value
maxN2 <- function(x, N = 2){
  len <- length(x)
  if(N > len){
    warning('N greater than length(x).  Setting N = length(x)')
    N <- length(x)
  }
  sort(x,partial = len - N+1)[len - N+1]
}


#Load covearge table
HPV16_coverage_path <- c("HPV16.coverage.csv")
coverage <- read.csv(paste(HPV16_coverage_path), sep="\t")

#counts coverage table
HPV16_counts_coverage_path<-c("/HPV16_long_counts_coverage_combined.csv")
counts_coverage <- read.csv(paste(HPV16_counts_coverage_path), header = T, sep = ",") 
counts_coverage$sample_id2<- counts_coverage$sample_id
coverage$sample_id2<- coverage$sample_id
coverage$sample_id2<-as.factor(coverage$sample_id2)

#fix coverage table
coverage$sample_id2 <- coverage$sample_id %>% {gsub(".*_[0-9]*-","", .)} %>% {gsub("_.*","", .)} %>% {gsub("-F_*$|-R_*$","",.)}

#filter out samples not included in the analysis
coverage<- coverage[!grepl("CIN|VIN|Pos|pos|Neg|neg|H2O", coverage$sample_id2),]
counts_coverage<- counts_coverage[!grepl("CIN|VIN|Pos|pos|Neg|neg|H2O", counts_coverage$sample_id2),]
coverage <- coverage[(coverage$sample_id2 %in% counts_coverage$sample_id2),]



# Filter out bases with QV <x
cutoff_qual<-30
coverage$A[coverage$qA < cutoff_qual] <- 0
coverage$T[coverage$qT < cutoff_qual] <- 0
coverage$C[coverage$qC < cutoff_qual] <- 0
coverage$G[coverage$qG < cutoff_qual] <- 0

# Filter out bases seen less than x times
coverage$A[coverage$A < 3] <- 0
coverage$T[coverage$T < 3] <- 0
coverage$C[coverage$C < 3] <- 0
coverage$G[coverage$G < 3] <- 0

# Separate F/R
Fvariants <- coverage[grep("-F", coverage$sample_id),]
Rvariants <- coverage[grep("-R", coverage$sample_id),]
Fvariants$sample_id<-gsub("-F","",Fvariants$sample_id)
Rvariants$sample_id<-gsub("-R","",Rvariants$sample_id)

#combine F and R
all_coverage=rbind(Fvariants, Rvariants) %>%
  group_by(sample_id2,chr,position, reference) %>% mutate(CumCoverage=sum(coverage),
                                                         cumA = sum(A),
                                                         cumG =sum(G),
                                                         cumC = sum(C),
                                                         cumT = sum (T)) %>%
  dplyr::select(sample_id2, chr, position, cumA, cumG, cumC, cumT,CumCoverage) %>%
  distinct()
#rename colnames
colnames(all_coverage) <- c("reference", "sample_id2", "chr", "position", "A", "G", "C", "T", "CumCoverage")


write.csv(all_coverage, file="/all_coverage_HPV16.csv")

# call variants with functions
all_coverage[, "major_base"] <- apply(all_coverage[5:8], 1, function (x) maxN(x, reachback = 0))
all_coverage[, "major_base_count"] <- apply(all_coverage[, 5:8],1, max)
all_coverage[, "minor_base"] <- apply(all_coverage[5:8], 1, function (x) maxN(x, reachback = 1))
all_coverage[, "minor_base_count"] <- apply(all_coverage[, 5:8],1, function (x) maxN2(x))
all_coverage$coverage <- apply((all_coverage[5:8]), 1 ,sum)
all_coverage$minor_base_freq <- all_coverage$minor_base_count / all_coverage$coverage *100
all_coverage$minor_base_freq[is.nan(all_coverage$minor_base_freq)] <- 0
all_coverage$minor_base[all_coverage$minor_base_count == 0] <- "N"

all_coverage$reference <- toupper(all_coverage$reference)



#applying cutoffs for iSNV calling
all_variants<-all_coverage %>% filter((coverage>100 & coverage <= 500 & minor_base_freq >= 5) | 
                                        (coverage>500 & coverage <= 1500 & minor_base_freq >= 1) |(coverage>1500 & coverage <= 10000 & minor_base_freq >= 1) |
                                        (coverage>10000 & minor_base_freq >= 1))
all_variants$sample_id2<-as.factor(all_variants$sample_id2)


# Remove homopolymeric regions within the non-coding and upstream regulatory regions

all_variants <- filter (all_variants, !(position >= 7282 & position <= 7286 |
                          position >= 7398 & position <= 7402 |
                          position >= 7472 & position <= 7476 |
                          position >= 7490 & position <= 7496 |
                          position >= 7605 & position <= 7610 |
                          position >= 7711 & position <= 7714 |
                          position >= 4158 & position <= 4175 |
                          position >= 4185 & position <= 4214))




write.csv(all_variants, file="/all_variants_HPV16.csv")

###### Pi calculated only among the variables positions (for test only and not used in the manuscript)

all_variants<-all_variants %>% mutate(depth = major_base_count + minor_base_count) %>% 
  mutate(pi_diff = major_base_count*minor_base_count) %>% mutate(pi_n_comp = (depth*(depth-1))/2) 
all_variants$pi<-all_variants$pi_diff/all_variants$pi_n_comp

##### pi per samples
pi_per_sample<-all_variants %>% group_by(sample_id2) %>% summarise(mean_pi = mean(pi)) 
var_per_sample<-all_variants %>% group_by(sample_id2) %>% tally()

pi_n_mnv_per_sample<- merge(x = pi_per_sample, y = var_per_sample,by.x = "sample_id2", by.y = "sample_id2", all.x = TRUE)

write.csv(pi_n_mnv_per_sample, file="/pi_n_isnv_per_sample_HPV16.csv")


#count the number of MNs per sample

var_per_sample<-all_variants %>% group_by(sample_id2) %>% tally() %>% rename(nr_MNVs = n, Sample_ID = sample_id2)

#### adding the trinucleotid context to the MNV_table from each sample separately
all_coverage$major_base[all_coverage$major_base_count < 20] <- NA

#transfer to anther table
all_coverage_major<-all_coverage

all_major<-all_coverage_major[,c("sample_id2", "position", "reference", "major_base")]
all_major$sample_id<-NULL
all_major$chr<-NULL

all_major$sample_id2<-as.factor(all_major$sample_id2)


s<-unique(all_major$sample_id2)
for (i in s) {
  assign(paste("s", i, sep = "_"), all_major[all_major$sample_id2 == i,] %>% 
           mutate(major_called = ifelse(is.na(major_base), paste(reference), paste(major_base)))) 
  
}


env <- .GlobalEnv

Names <- ls(pattern = "^s_[0-9]*-", env)

L <- mget(Names, env)

all_major_bound<-do.call("rbind", L)

rm(list=ls(pattern="^s_[0-9]*-"))



all_major_bound<-all_major_bound %>% mutate(next_base = c(""), previous_base = c(""))%>% group_by(sample_id2, position)


for (i in unique(all_major_bound$sample_id2)) {
  all_major_bound$previous_base[all_major_bound$sample_id2 == i]<-lag(all_major_bound$major_called[all_major_bound$sample_id2 ==i])
  all_major_bound$next_base[all_major_bound$sample_id2 == i]<-lead(all_major_bound$major_called[all_major_bound$sample_id2 ==i])
  all_major_bound[all_major_bound$sample_id2 == i,]<-all_major_bound[all_major_bound$sample_id2 == i,] %>% mutate_at(vars(previous_base), ~replace_na(.,tail(major_called, n=1))) %>%
    mutate_at(vars(next_base), ~replace_na(.,head(major_called, n=1)))
}

all_major_bound<-all_major_bound %>% select(sample_id2,position, previous_base,major_called,next_base)
all_major_bound$major_base <- all_major_bound$major_called
all_major_bound$major_called<-NULL
major_selected_samples<-unique(all_major_bound$sample_id2)
all_variants<-all_variants[all_variants$sample_id2 %in% major_selected_samples,]

variants_called_merged<-left_join(all_variants, all_major_bound, by = c("sample_id2","major_base", "position"))


#adding mutatiotonal signature
variants_called_merged$mutation<-paste(variants_called_merged$major_base, ">", variants_called_merged$minor_base, sep = "")


#converting mutational signature into 6 categories 

variants_called_merged$com_previous_base <- ifelse(variants_called_merged$previous_base == "A", "T",
                                                   ifelse(variants_called_merged$previous_base == "T", "A",
                                                          ifelse(variants_called_merged$previous_base == "G", "C",
                                                                 "G")))
variants_called_merged$com_next_base <- ifelse(variants_called_merged$next_base == "A", "T",
                                               ifelse(variants_called_merged$next_base == "T", "A",
                                                      ifelse(variants_called_merged$next_base == "G", "C",
                                                             "G")))
#make six categories instead of 12 for mutation

variants_called_merged$six_mut <-ifelse(variants_called_merged$mutation == "A>T", "T>A_1",
                                        ifelse(variants_called_merged$mutation == "A>C", "T>G_1",
                                               ifelse(variants_called_merged$mutation == "A>G", "T>C_1",
                                                      ifelse(variants_called_merged$mutation == "G>A", "C>T_1",
                                                             ifelse(variants_called_merged$mutation == "G>T", "C>A_1",
                                                                    ifelse(variants_called_merged$mutation == "G>C", "C>G_1", variants_called_merged$mutation))))))

variants_called_merged$tri_context<-ifelse(grepl("_1",variants_called_merged$six_mut), paste(variants_called_merged$com_previous_base, "*", variants_called_merged$com_next_base, sep = ""),
                                           paste(variants_called_merged$previous_base, "*", variants_called_merged$next_base, sep = ""))

variants_called_merged$six_mut<-gsub("_1", "", variants_called_merged$six_mut)

##Add genes
HPV_genes<-read.csv("./Clade9_subtypes_all.csv", header = TRUE)
HPV_genes<-HPV_genes[HPV_genes$SUB == HPV_subtype,]

variants_called_merged$gene<-ifelse(variants_called_merged$position >= HPV_genes$START[HPV_genes$GENE == "URR1"] & variants_called_merged$position <= HPV_genes$END[HPV_genes$GENE == "URR1"],"URR",
                                    ifelse(variants_called_merged$position >= HPV_genes$START[HPV_genes$GENE == "E6"] & variants_called_merged$position <= HPV_genes$END[HPV_genes$GENE == "E6"],"E6",
                                           ifelse(variants_called_merged$position >= HPV_genes$START[HPV_genes$GENE == "E7"] & variants_called_merged$position <= HPV_genes$END[HPV_genes$GENE == "E7"],"E7",
                                                  ifelse(variants_called_merged$position >= HPV_genes$START[HPV_genes$GENE == "E1"] & variants_called_merged$position <= HPV_genes$END[HPV_genes$GENE == "E1"],"E1",
                                                         ifelse(variants_called_merged$position >= HPV_genes$START[HPV_genes$GENE == "E2"] & variants_called_merged$position <= HPV_genes$END[HPV_genes$GENE == "E2"],"E2",
                                                                ifelse(variants_called_merged$position >= HPV_genes$START[HPV_genes$GENE == "E5"] & variants_called_merged$position <= HPV_genes$END[HPV_genes$GENE == "E5"],"E5",
                                                                       ifelse(variants_called_merged$position >= HPV_genes$START[HPV_genes$GENE == "NCR"] & variants_called_merged$position <= HPV_genes$END[HPV_genes$GENE == "NCR"],"NCR",
                                                                              ifelse(variants_called_merged$position >= HPV_genes$START[HPV_genes$GENE == "L2"] & variants_called_merged$position <= HPV_genes$END[HPV_genes$GENE == "L2"],"L2",
                                                                                     ifelse(variants_called_merged$position >= HPV_genes$START[HPV_genes$GENE == "L1"] & variants_called_merged$position <= HPV_genes$END[HPV_genes$GENE == "L1"],"L1",
                                                                                            ifelse(variants_called_merged$position >= HPV_genes$START[HPV_genes$GENE == "URR"] & variants_called_merged$position <= HPV_genes$END[HPV_genes$GENE == "URR"], "URR", "NA"))))))))))



variants_called_merged[,c("sample_id2", "mutation", "gene")] <- lapply(variants_called_merged[,c("sample_id2", "mutation", "gene")], factor)

length(unique(variants_called_merged$sample_id2))
tapply(variants_called_merged$mutation, variants_called_merged$sample_id2, summary)

write.csv(variants_called_merged, file="/variants_called_merged_HPV16.csv")



## Selection loci

all_variants_unique<-unique(all_variants$position)
all_coverage_select<-all_coverage %>% filter(position %in% all_variants_unique)

###### Pi calculated over the entire genome (used in the manuscript)

all_coverage_select<-all_coverage_select %>% mutate(depth = major_base_count + minor_base_count) %>% 
mutate(pi_diff = major_base_count*minor_base_count) %>% mutate(pi_n_comp = (depth*(depth-1))/2) 
all_coverage_select$pi<-all_coverage_select$pi_diff/all_coverage_select$pi_n_comp
all_coverage_select$pi[is.na(all_coverage_select$pi)]<-0

##### pi per samples
pi_per_sample_select<-all_coverage_select %>% group_by(sample_id2) %>% summarise(mean_pi = mean(pi)) 
var_per_sample_select<-all_coverage_select %>% group_by(sample_id2) %>% tally()

pi_n_mnv_per_sample_select<- merge(x = pi_per_sample_select, y = var_per_sample_select,by.x = "sample_id2", by.y = "sample_id2", all.x = TRUE)

write.csv(pi_n_mnv_per_sample_select, file="/pi_correct_HPV16.csv")



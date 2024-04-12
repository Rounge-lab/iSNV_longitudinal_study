rm(list=ls())
par(mfrow=c(1,1))

setwd("")
source("HighstatLibV6.R")
# detach("package:ape", unload = TRUE)

library("ggplot2")
library("ggpubr")
library("ggeffects")
library("purrr")
library("tidyr")
library("stringr")
library("reshape2")
library("insect")
library("seqinr")
library("adegenet")
library("tools")
library("latexpdf")
library("viridis")
library("MASS")
library("gridExtra")
library("multcomp")
library("lattice")
library("car")
library("Matrix")
library("glmmTMB")
library("DHARMa")
library("compare")
library("dplyr")
library("bbmle")
library("geepack")
library("lme4")
library("lmerTest")
library("aods3")
library("performance")
library("rstatix")


##### Data setup

# Data from MNV_calling_long.R

pi_n_mnv_per_sample <- read.csv("./pi_n_mnv_per_sample_06_2023.csv")
pi_n_mnv_per_sample_v2 <- read.csv("./pi_n_mnv_per_sample_V2_06_2023.csv")
freq_per_sample <-read.csv("./freq_per_sample.csv") 
pi_n_mnv_per_sample$sample_id2<-as.factor(pi_n_mnv_per_sample$sample_id2)
pi_n_mnv_per_sample_v2$sample_id2<-as.factor(pi_n_mnv_per_sample_v2$sample_id2)
variants_called<- read.csv("./variants_called_merged_06_2023.csv")


# Data from script: Merge_count_coverage_diagnosis.R
df_total <- read.csv("./Merged_cov_diag_total_ME_Redit_10_2023.csv")



# Merge datasets
df_total2<-merge(x=df_total, y=pi_n_mnv_per_sample, by.x="sample_id2", by.y="sample_id2", all.x=TRUE)
df_total2.5<-merge(x=df_total2, y=pi_n_mnv_per_sample_v2, by.x="sample_id2", by.y="sample_id2", all.x=TRUE)
df_total3<-merge(x=df_total2.5, y=freq_per_sample, by.x="sample_id2", by.y="sample_id2", all.x=TRUE)

variants_called_merged<-merge(x=variants_called, y=df_total, by.x="sample_id2", by.y="sample_id2", all.x=TRUE)


## Remove samples that have been found with coinfection
df_total3 <- filter(df_total3, sample_id2 != "sample1-HPV16" & sample_id2 != "sample2-HPV16" & sample_id2 != "sample3-HPV16" & sample_id2 != "sample4-HPV16")
df_total3$time_diff_fromA2<-as.integer(df_total3$time_diff_fromA2)

df_total3_order <- df_total3[order(df_total3$Personnummer,df_total3$time_diff_fromA2),]


################################# Data summary ###############################

summary_diag <- df_total3_order %>%
  group_by(Grouping_2023) %>%
  summarize(
    merge2 = n(),
    Personnummer = n_distinct(Personnummer)
  )

summary_total <- df_total3_order %>%
  summarize(
    Grouping_2023 = "Total",
    merge2 = n(),
    Personnummer = n_distinct(Personnummer)
  )

summary_dataset <- bind_rows(summary_diag, summary_total)
summary_dataset


### Total summary mean coverage and  percent genome covered 100x 
summary(df_total3_order[c("mean_coverage", "genome_covered_100x")])





#####################################################################
############################### Model number iSNVs ##################



# fit the model

glmm_mnv <- glmmTMB(n.x ~  Grouping_2023 * time_diff_fromA2 + (1 | Personnummer), family = nbinom2, data = df_total3_order)


# print the summary
summary(glmm_mnv)

# check the confidence interval
confint(glmm_mnv)


## Investigate residuals with DHARMA

# simulationOutput <- simulateResiduals(fittedModel = glmm_mnv, plot=T)
simulationOutput <- simulateResiduals(fittedModel = glmm_mnv, n=1000, plot=T)

testResiduals(simulationOutput)

# Test overdispersion
testDispersion(simulationOutput)

# Test zeroinflation
testZeroInflation(simulationOutput)


# Plot the predicted relationship

mydf_mnv <- ggpredict(glmm_mnv, terms = c("time_diff_fromA2", "Grouping_2023"),  ci.lvl = 0.95)
mydf_mnv$group<-factor(mydf_mnv$group, levels =c("High_grade", "Low_grade", "Normal"))
mydf_mnv_order<- mydf_mnv %>% arrange(x, group)
mydf_mnv_subset<- mydf_mnv_order %>% filter((group=="High_grade") & x <= 700 | (group=="Low_grade") & x <= 700 | (group=="Normal"))

ggplot(mydf_mnv_subset, aes(x, predicted)) +
  scale_color_manual(values=c("#440154", "#21918c", "#fde725"), guide = "none")+
  scale_fill_manual(values=c( "#440154", "#21918c", "#fde725"), name = "Diagnostic category", labels = c("High grade", "Low grade", "Normal" ))+
  theme_bw() +
  theme(axis.text=element_text(size=18),
        axis.title=element_text(size=16,face="bold"),
        legend.title = element_text(size=18),
        legend.text = element_text(size=16))+
  theme(panel.grid.major = element_blank())+
  geom_line(aes(color = group), size=1.5)+
  geom_ribbon(data=mydf_mnv_subset,aes(ymin = conf.low, ymax = conf.high, fill = group), alpha = 0.08, color="grey", size=0.1) +
  labs(y = "Number of iSNV", x = "Time difference (days)")+
  geom_point(data=df_total3_order, aes(time_diff_fromA2, n.x, fill = Grouping_2023), alpha = 0.8, size=4.2,pch=21)+
  guides(fill = guide_legend(override.aes = list(size=10)))


########## Generate summary iSNV #############

df_summary_iSNV <- df_total3_order %>%
  group_by(Grouping_2023) %>%
  summarise(
    Mean_iSNV= mean(n.x),
    Min_iSNV = min(n.x),
    Max_iSNV = max(n.x))




#######################################################################
###################### Model mean pi ##################################

lmer_pi<- lmer((mean_pi.y) ~ Grouping_2023 * time_diff_fromA2 + (1 | Personnummer), data = df_total3_order)

# print the summary
summary(lmer_pi)

# check the confidence interval
confint(lmer_pi)

## Investigate residuals with DHARMA

# simulationOutput <- simulateResiduals(fittedModel = glmm_mnv, plot=T)
simulationOutput <- simulateResiduals(fittedModel = lmer_pi, n=1000, plot=T)

testResiduals(simulationOutput)

# Test overdispersion
testDispersion(simulationOutput)

# Test zeroinflation
testZeroInflation(simulationOutput)


# Plot the predicted relationship
mydf_pi <- ggpredict(lmer_pi, terms = c("time_diff_fromA2", "Grouping_2023"),  ci.lvl = 0.95)
mydf_pi$group<-factor(mydf_pi$group, levels =c("High_grade", "Low_grade", "Normal"))
mydf_pi_order<- mydf_pi %>% arrange(x, group)
mydf_pi_subset<- mydf_pi_order %>% filter((group=="High_grade") & x <= 400 | (group=="Low_grade") & x <= 600 | (group=="Normal"))


ggplot(mydf_pi_subset, aes(x, predicted)) +
  scale_color_manual(values=c("#440154", "#21918c", "#fde725"), guide = "none")+
  scale_fill_manual(values=c("#440154", "#21918c", "#fde725"), name = "Diagnostic category", labels = c("High grade", "Low grade", "Normal"))+
  theme_bw() +
  theme(axis.text=element_text(size=18),
        axis.title=element_text(size=16,face="bold"),
        legend.title = element_text(size=18),
        legend.text = element_text(size=16))+
  theme(panel.grid.major = element_blank())+
  geom_line(aes(color = group), size=1.5)+
  scale_y_continuous(limits = c(0.0005,0.004))+
  geom_ribbon(data=mydf_pi_subset,aes(ymin = conf.low, ymax = conf.high, fill = group), alpha = 0.08, color="grey", size=0.1) +
  labs(y = "Nucleotide diversity (\u03C0)", x = "Time difference (days)")+
  geom_point(data=df_total3_order, aes(time_diff_fromA2, mean_pi.y, fill = Grouping_2023), alpha = 0.8, size=4.2,pch=21)+
  guides(fill = guide_legend(override.aes = list(size=10)))


########## Generate summary data quality table #############

df_summary_pi <- df_total3_order %>%
  group_by(Grouping_2023) %>%
  summarise(
    Mean_pi= mean(mean_pi.y),
    Min_pi = min(mean_pi.y),
    Max_pi = max(mean_pi.y),
    Var_pi = var(mean_pi.y))

df_summary_pi




###################################################################################
###################### Model mutation signature ##################################


##################### Mutation signature Graph




signature_all2 <- variants_called_merged %>%
  group_by(sample_id2, Grouping_2023, six_mut, tri_context) %>%
  summarise(count = n())


##################### Mutation signature models


## Preparing the dataset
signature_ct_TCA<-subset(signature_all2, six_mut=="C>T" & tri_context=="T*A") 
signature_ct_TCA <- signature_ct_TCA[,-2]

signature_ct_TCA_TCT<-subset(signature_all2, six_mut=="C>T" & tri_context=="T*A" | six_mut=="C>T" & tri_context=="T*T") 
signature_ct_TCA_TCT2 <- signature_ct_TCA_TCT%>%
  group_by(sample_id2, six_mut) %>%
  summarise(count_TCA_TCT = sum(count)) 



df_total3_signature<-merge(x=df_total3_order, y=signature_ct_TCA , by.x="sample_id2", by.y="sample_id2", all.x=TRUE)
df_total3_signature<-merge(x=df_total3_signature, y=signature_ct_TCA_TCT2 , by.x="sample_id2", by.y="sample_id2", all.x=TRUE)
df_total3_signature$count[is.na(df_total3_signature$count)] <- 0
df_total3_signature$count_TCA_TCT[is.na(df_total3_signature$count_TCA_TCT)] <- 0
df_total3_order_signature<- df_total3_signature[order(df_total3_signature$Personnummer,df_total3_signature$time_diff_fromA2),]



sum(df_total3_order_signature$count_TCA_TCT == 0)



# Corrected for zeroinflation
glmmtmb_model_zi<- glmmTMB(count_TCA_TCT ~ Grouping_2023 * time_diff_fromA2 + (1 | Personnummer), family = poisson, ziformula=~1, data = df_total3_order_signature)



# print the summary
summary(glmmtmb_model_zi)


# check the confidence interval
confint(glmmtmb_model_zi)

## Investigate residuals with DHARMA

# simulationOutput <- simulateResiduals(fittedModel = glmm_mnv, plot=T)
simulationOutput <- simulateResiduals(fittedModel = glmmtmb_model_zi, n=1000, plot=T)

testResiduals(simulationOutput)

# Test overdispersion
testDispersion(simulationOutput)

# Test zeroinflation
testZeroInflation(simulationOutput)


## Prepare plot

mydf_mutsig <- ggpredict(glmmtmb_model_zi, terms = c("time_diff_fromA2", "Grouping_2023"),  ci.lvl = 0.95)
mydf_mutsig$group<-factor(mydf_mutsig$group, levels =c("High_grade", "Low_grade", "Normal"))
mydf_mutsig_order<- mydf_mutsig %>% arrange(x, group)
mydf_mutsig_subset<- mydf_mutsig_order %>% filter((group=="High_grade") & x <= 700 | (group=="Low_grade") & x <= 700 | (group=="Normal"))


ggplot(mydf_mutsig_subset, aes(x, predicted)) +
  scale_color_manual(values=c("#440154", "#21918c", "#fde725"), guide = "none")+
  scale_fill_manual(values=c("#440154", "#21918c", "#fde725"), name = "Diagnostic category", labels = c("High grade", "Low grade", "Normal"))+
  theme_bw() +
  theme(axis.text=element_text(size=18),
        axis.title=element_text(size=16,face="bold"),
        legend.title = element_text(size=18),
        legend.text = element_text(size=16))+
  theme(panel.grid.major = element_blank())+
  geom_line(aes(color = group), size=1.5)+
  geom_ribbon(data=mydf_mutsig_subset,aes(ymin = conf.low, ymax = conf.high, fill = group), alpha = 0.08, color="grey", size=0.1) +
  labs(y = "Number of C>T mutations \nin T*A and T*T context", x = "Time difference (days)")+
  geom_point(data=df_total3_order_signature, aes(time_diff_fromA2, count_TCA_TCT, fill = Grouping_2023), alpha = 0.8, size=4.2,pch=21)+
  guides(fill = guide_legend(override.aes = list(size=10)))


########## Generate summary data quality table #############

df_summary_APOBEC <- df_total3_order_signature %>%
  group_by(Grouping_2023) %>%
  summarise(
    Mean_n_mut= mean(count_TCA_TCT),
    Min_n_mut = min(count_TCA_TCT),
    Max_n_mut = max(count_TCA_TCT))




########################################################
########## Models without accountig for time ###########
########################################################



################ Model number iSNVs ############### 


# fit the model
n.x <- df_total3_order$n.x
Grouping_2023 <- as.factor(df_total3_order$Grouping_2023)
Personnummer <- as.factor(df_total3_order$Personnummer)

glmm_mnv_no_time <- glmmTMB(n.x ~  Grouping_2023 + (1 | Personnummer), family = nbinom2 )


# print the summary
summary(glmm_mnv_no_time)

# check the confidence interval
confint(glmm_mnv_no_time)

# Check residuals
simulationOutput <- simulateResiduals(fittedModel = glmm_mnv_no_time, n=1000, plot=T)
testResiduals(simulationOutput)
# Test overdispersion
testDispersion(simulationOutput)
# Test zeroinflation
testZeroInflation(simulationOutput)



    
  


###################### Model mean pi ####################


lmer_pi_no_time<- lmer(mean_pi.y ~ Grouping_2023+ (1 | Personnummer), data = df_total3_order)

# print the summary
summary(lmer_pi_no_time)

# check the confidence interval
confint(lmer_pi_no_time)

## Investigate residuals with DHARMA

# simulationOutput <- simulateResiduals(fittedModel = glmm_mnv, plot=T)
simulationOutput <- simulateResiduals(fittedModel = lmer_pi_no_time, n=1000, plot=T)

testResiduals(simulationOutput)

# Test overdispersion
testDispersion(simulationOutput)

# Test zeroinflation
testZeroInflation(simulationOutput)





############### Model mutation signature ##############


glmmtmb_model_zi_no_time<- glmmTMB(count_TCA_TCT ~ Grouping_2023 + (1 | Personnummer), family = poisson, ziformula=~1, data = df_total3_order_signature)



# print the summary
summary(glmmtmb_model_zi_no_time)

# check the confidence interval
confint(glmmtmb_model_zi_no_time)

## Investigate residuals with DHARMA

# simulationOutput <- simulateResiduals(fittedModel = glmm_mnv, plot=T)
simulationOutput <- simulateResiduals(fittedModel = glmmtmb_model_zi_no_time, n=1000, plot=T)

testResiduals(simulationOutput)

# Test overdispersion
testDispersion(simulationOutput)

# Test zeroinflation
testZeroInflation(simulationOutput)





########################################################
################# Supplementary plots ##################
########################################################


# Relationship between number of iSNV and diagnosis grouping
plot_grid<-ggplot(df_total3_order, aes(x=time_diff_fromA2, y=n.x)) +
  theme_bw()+
  # geom_smooth(method='lm', aes(color = NA), formula=y~0+x, se=TRUE) +
  geom_point(aes(color=Grouping_2023), size=6) + geom_line() +
  scale_color_manual(values=c("#440154", "#21918c",  "#fde725"), name = "Diagnostic category", labels = c("High grade", "Low grade", "Normal" ))+
  # scale_color_manual(values=c( "#fde725", "#21918c",  "#440154"), name = "Diagnostic category", labels = c("High grade", "Low grade", "Normal" ))+
  facet_wrap(~ Personnummer, ncol = 7) +
  theme(strip.background = element_blank(),
        strip.text.x.top  = element_blank())+
  labs(x="Time difference (days)", y='Number of iSNV') +
  theme(axis.text=element_text(size=14),
        axis.title=element_text(size=16,face="bold"),
        legend.title = element_text(size=16),
        legend.text = element_text(size=14))
plot_grid


# Relationship between Pi and diagnosis grouping
plot_grid<-ggplot(df_total3_order, aes(x=time_diff_fromA2, y=mean_pi.y)) +
  theme_bw()+
  # geom_smooth(method='lm', aes(color = NA), formula=y~0+x, se=TRUE) +
  geom_point(aes(color=Grouping_2023), size=6) + geom_line() +
  scale_color_manual(values=c("#440154", "#21918c",  "#fde725"), name = "Diagnostic category", labels = c("High grade", "Low grade", "Normal" ))+
  # scale_color_manual(values=c( "#fde725", "#21918c",  "#440154"), name = "Diagnostic category", labels = c("High grade", "Low grade", "Normal" ))+
  facet_wrap(~ Personnummer, ncol = 7) +
  theme(strip.background = element_blank(),
        strip.text.x.top  = element_blank())+
  labs(x="Time difference (days)", y="Nucleotide diversity (\u03C0)") +
  theme(axis.text=element_text(size=14),
        axis.title=element_text(size=16,face="bold"),
        legend.title = element_text(size=16),
        legend.text = element_text(size=14))
plot_grid


# Relationship between APOBEC3 and diagnosis grouping
plot_grid<-ggplot(df_total3_order_signature, aes(x=time_diff_fromA2, y=count_TCA_TCT)) +
  theme_bw()+
  # geom_smooth(method='lm', aes(color = NA), formula=y~0+x, se=TRUE) +
  geom_point(aes(color=Grouping_2023), size=6) + geom_line() +
  scale_color_manual(values=c("#440154", "#21918c",  "#fde725"), name = "Diagnostic category", labels = c("High grade", "Low grade", "Normal" ))+
  # scale_color_manual(values=c( "#fde725", "#21918c",  "#440154"), name = "Diagnostic category", labels = c("High grade", "Low grade", "Normal" ))+
  facet_wrap(~ Personnummer, ncol = 7) +
  theme(strip.background = element_blank(),
        strip.text.x.top  = element_blank())+
  labs(x="Time difference (days)", y="Number of C>T mutations \nin T*A and T*T context") +
  theme(axis.text=element_text(size=14),
        axis.title=element_text(size=16,face="bold"),
        legend.title = element_text(size=16),
        legend.text = element_text(size=14))
plot_grid


#################


## Review test model with age



glmmtmb_model_iSNV_age<- glmmTMB(n.x ~ alder + (1 | Personnummer), family =nbinom2, data = df_total3_order_signature)
glmmtmb_model_iSNV_age2<- glmmTMB(n.x ~ Grouping_2023 * alder + (1 | Personnummer), family =nbinom2, data = df_total3_order_signature)

glmmtmb_model_APOBEC_age<- glmmTMB(count_TCA_TCT ~ alder + (1 | Personnummer), family = poisson, ziformula=~1, data = df_total3_order_signature)
glmmtmb_model_APOBEC_age2<- glmmTMB(count_TCA_TCT ~ Grouping_2023 * alder + (1 | Personnummer), family = poisson, ziformula=~1, data = df_total3_order_signature)

summary(glmmtmb_model_iSNV_age)
summary(glmmtmb_model_APOBEC_age)
summary(glmmtmb_model_iSNV_age2)
summary(glmmtmb_model_APOBEC_age2)





# Count the number of participant with the same or different diagnostic over time
diagnostic_counts <- df_total3_order %>%
  group_by(Personnummer) %>%
  summarize(num_diagnostic_categories = n_distinct(Grouping_2023))

# same diagnostic category
same_diagnostic_count <- diagnostic_counts %>%
  filter(num_diagnostic_categories == 1) %>%
  nrow()

# different diagnostic categories
different_diagnostic_count <- diagnostic_counts %>%
  filter(num_diagnostic_categories > 1) %>%
  nrow()

# Print the results
cat("Number of patients with all samples in the same diagnostic category:", same_diagnostic_count, "\n")
cat("Number of patients with samples in different diagnostic categories:", different_diagnostic_count, "\n")


# Create table with minor base frequency, mean pi and number if iSNV  for each diagnostic group
aggregate(df_total3_order$minor_base_freq~df_total3_order$Grouping_2023, FUN=mean)
aggregate(df_total3_order$mean_pi.y~df_total3_order$Grouping_2023,  FUN=mean)
aggregate(df_total3_order$n.x~df_total3_order$Grouping_2023,  FUN=mean)


#=============================================================================#
# Script created to import data and conduct analyses for
# LV Barbosa et al. Austral Ecology 
# Script compatible with R version R version 3.6.3 (2020-02-29)
#=============================================================================#

##############################################################

# Load required libraries
library(tidyverse)
library(lme4)
library(xtable)
##############################################################
# We used split-plot analysis of variance (ANOVA) to determine how N-fertilization 
# influenced litter decomposition. The main effects were 
# litter origin (i.e., High N, low N, or control plots) and 
# litter destination (i.e., High N, low N, or control plots);
# we also included the litter origin x litter destination interaction. 
# Data on the proportion of litter biomass remaining was logit-transformed prior to analysis
# Separate analyses were conducted for each year of the experiment. 

##############################################################
# LOAD AND PROCESS DATA: litter decomposition
##############################################################
decomp_data<-read.csv("./Data/LVB_decomp_data_1march2018.csv", dec=".", header = TRUE, sep = ",", check.names=FALSE )
cols <- c("key", "number", "plot")
decomp_data[cols]<-lapply(decomp_data[cols],as.factor)
rm(cols)
decomp_data <- decomp_data %>% mutate(destination = ordered(destination, levels = c("Control","Low-N","High-N")))
decomp_data <- decomp_data %>% mutate(source = ordered(source, levels = c("Control","Low-N","High-N")))
# calculate the proportion lost from proportion remaining
decomp_data$prop.lost<-(1-decomp_data$prop.rem.3)
hist(decomp_data$prop.lost)
#logit transformation of prop lost 
decomp_data$logit_decomp = log10(decomp_data$prop.lost/(1 - decomp_data$prop.lost))
str(decomp_data)
summary(decomp_data)
decomp_summary_dest <- decomp_data %>% 
  group_by(source,destination,year,replicate, replicate2)
decomp_summary_dest
str(decomp_summary_dest)

# TABLE OF OVERALL DECOMPOSITION BY YEAR (MEAN +/- SD)
decomp.yr <-decomp_data %>% select(year,prop.lost) %>% 
  group_by(year) %>%  
  summarize(mean = mean(prop.lost, na.rm = T), SD = sd(prop.lost, na.rm = T))

# Need rto calclulate averages of logit decomp in each plot for split-plot anova
# first calculate the average SD and SE of the mean of the averages of the 5 bags from each treatement that have been placed in each plot
plot.data<-ungroup(decomp_summary_dest)  
plot.data_M1<-plot.data %>% select(year, source, destination, plot, 
                                   prop.rem.3, prop.lost,logit_decomp) %>%
  group_by(year,source,destination,plot) %>% 
  summarize(avg.prop.lost = mean(prop.lost, na.rm = T),
            avg.logit.decomp = mean(logit_decomp, na.rm = T)) 

plot.data_M2<-plot.data_M1 %>%  group_by(year,source,destination) %>% 
  summarize(avg.prop.lost = mean(avg.prop.lost, na.rm = T),avg.logit.decomp = mean(avg.logit.decomp, na.rm = T))

plot.data_M3<-plot.data_M1 %>%  group_by(year,source,destination) %>% 
  summarize(SD.prop.lost = sd(avg.prop.lost, na.rm = T),SE.prop.lost = sd(avg.prop.lost, na.rm = T)/sqrt(length(avg.prop.lost)),SD.logit.decomp = sd(avg.logit.decomp, na.rm = T), SE.logit.decomp = sd(avg.logit.decomp, na.rm = T)/sqrt(length(avg.logit.decomp)))

plot.data_M4<-left_join(plot.data_M2,plot.data_M3,by=c("year","source","destination"))

plot.data<-left_join(plot.data,plot.data_M4,by=c("year","source","destination"))
avg.2009<-filter(plot.data_M1,year==2009)
avg.2010<-filter(plot.data_M1,year==2010)
avg.2011<-filter(plot.data_M1,year==2011)

DATA.AVG<-avg.2009
# DATA.AVG<-avg.2010
# DATA.AVG<-avg.2011

# SPLIT PLOT ANOVA http://www.personal.psu.edu/mar36/stat_461/split_plot/split_plot.html

# 2009
sp.aov.2009 <- aov(avg.logit.decomp ~ source*destination+Error(plot),data = avg.2009)
summary(sp.aov.2009)

# 2010
sp.aov.2010 <- aov(avg.logit.decomp ~ source*destination+Error(plot),data = avg.2010)
summary(sp.aov.2010)

# 2011
sp.aov.2011 <- aov(avg.logit.decomp ~ source*destination+Error(plot),data = avg.2011)
summary(sp.aov.2011)

# OUTPUT OF AOV TABLES TO WORD DOC
# 2009
out.2009<-xtable(summary(sp.aov.2009))
out.2009<-round(out.2009, 4)
write.csv(out.2009, file="./Output/sp-aov-2009.csv",na = "")
# 2010
out.2010<-xtable(summary(sp.aov.2010))
out.2010<-round(out.2010, 4)
write.csv(out.2010, file="./Output/sp-aov-2010.csv",na = "")
# 2011
out.2011<-xtable(summary(sp.aov.2011))
out.2011<-round(out.2011, 4)
write.csv(out.2011, file="./Output/sp-aov-2011.csv",na = "")

# open word, import file to document
# select the txt file
# once imported, highlight and under TABLE tab, convert text to table with comma sep.



##############################
# FIGURE 
##############################

plot.data_M4$destination <- factor(plot.data_M4$destination, levels = c("High-N","Low-N","Control"))
plot.data_M4$source <- factor(plot.data_M4$source, levels = c("High-N","Low-N","Control"))
fig2 <- ggplot(plot.data_M4,aes(x=source, y=avg.prop.lost, fill=destination)) +
  facet_grid(. ~ year)+
  geom_bar(stat="identity", position=position_dodge(), color="black") +
  geom_errorbar(aes(ymin=avg.prop.lost-SE.prop.lost , ymax=avg.prop.lost+SE.prop.lost ), width=.3,
                position=position_dodge(.9))+
  scale_y_continuous(limit=c(0.0, 0.5))+
  xlab("Litter Source")+
  ylab("Mass loss (proportion)")+
  labs(fill=guide_legend(title="Litter Destination"))+
  scale_fill_manual(values = c("Control" = "#FFFFFF", "Low-N" = "#CCCCCC", "High-N" = "#666666"))
#   # scale_colour_manual(values=c("#666666","#0072B2","#000066"))
fig2<-fig2 + theme_classic() + theme(panel.border = element_blank(), panel.grid.major = element_blank(), 
                                     axis.line.y = element_line(color="black", size = 0.5, lineend="square"),
                                     axis.line.x = element_line(color="black", size = 0.5, lineend="square"),
                                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), #sets colors of axes
                                     plot.title = element_text(hjust=0.05, vjust=-1.8, face="bold", size=22),        #Sets title size, style, location
                                     axis.title.x=element_text(colour="black", size = 20, vjust=-2),            #sets x axis title size, style, distance from axis #add , face = "bold" if you want bold
                                     axis.title.y=element_text(colour="black", size = 20, vjust=2),            #sets y axis title size, style, distance from axis #add , face = "bold" if you want bold
                                     legend.position = c(1.02, 0.75),
                                     # legend.title = element_blank(),
                                     strip.text = element_text(colour="black", size = 16,face="bold"),                    
                                     axis.text.x=element_text(colour="black", size = 16, angle=320,hjust=-0.1),                              #sets size and style of labels on axes
                                     axis.text.y=element_text(colour="black", size = 16),                              #sets size and style of labels on axes
                                     plot.margin = unit(c(1,3,2,1), "cm"),
                                     panel.spacing = unit(1, "lines"))

fig2

# 900 x 500 to save
##############################################################

# We used nested ANOVA to test how N-fertilization affect N litter concentration. 
# Levels of N addition where the litter was collected (litter origin) and the year of 
# litter harvested were fixed effects, with the two measurements of each composite nested 
# within the plot in which the litter was collected (note that we did not use repeated measures 
# ANOVA because we did not necessarily sample the same individual plant each year). 

##############################################################
# UPLOAD & STANDARDIZE DATA on Litter percent nitrogen
##############################################################
litter_N_data<-read.csv("./Data/litter_N.csv", dec=".", header = TRUE, sep = ",", check.names=FALSE )
litter_N_data$treatment<-recode(litter_N_data$treatment, 'high-N' = "High-N",'low-N' = "Low-N",'control'="Control")
litter_N_data$measurement<-rep(1:2, 45)
litter_N_data[8]<-NULL
litter_N_data$year<-as.factor(litter_N_data$year)
str(litter_N_data)
hist(litter_N_data$N_percent)
litter_N_data$log_N_percent<-log(litter_N_data$N_percent)
hist(litter_N_data$log_N_percent)

#This calclulates the mean percent N of the two samples from each plot in each year
litter_n_plot_avg<-litter_N_data %>% group_by(year,treatment,plot) %>% summarize(mean_percent_N=mean(N_percent))
litter_n_plot_avg<-litter_n_plot_avg %>% rename("source"="treatment")
#This calclulates the mean percent N of all samples from each trt in each year
litter_n_trt_avg<-litter_N_data %>% group_by(year,treatment) %>% summarize(mean_percent_N=mean(N_percent))
litter_n_trt_avg<-litter_n_trt_avg %>% rename("source"="treatment")
litter_n_trt_avg$year<-as.numeric(levels(litter_n_trt_avg$year))[litter_n_trt_avg$year]
litter_n_trt_avg <- litter_n_trt_avg %>% mutate(source = ordered(source, levels = c("Control","Low-N","High-N")))
litter_n_trt_avg <- as.data.frame(litter_n_trt_avg)

#############################
# analyses
#############################
# LEAF N BY TREATMENT
# have to avg the two measurments for each plot
lfN.trt <-litter_N_data %>% select(N_percent,plot,treatment,measurement) %>% 
  group_by(plot, treatment) %>% 
  summarize(mean = mean(N_percent, na.rm = T)) %>% 
  group_by(treatment) %>% 
  summarize(mean = mean(mean, na.rm = T))
lfN.trt

# have to avg the two measurments for each plot
lfN.yr <-litter_N_data %>% select(N_percent,plot,year,measurement) %>% 
  group_by(plot, year) %>% 
  summarize(mean = mean(N_percent, na.rm = T)) %>% 
  group_by(year) %>% 
  summarize(mean = mean(mean, na.rm = T))
lfN.yr

########################
#  ANOVA OF litter N
########################
litter_N_data_wide<-litter_N_data %>% spread(measurement, N_percent)
litter_N_data_avg <-litter_N_data %>% group_by(year, treatment,plot) %>% summarise(avg = mean(N_percent))
litter_N_data_avg$year<-as.factor(litter_N_data_avg$year)
litter_N_data$treatment <- factor(litter_N_data$treatment, levels = c("High-N","Low-N","Control"))

# https://stackoverflow.com/questions/37497948/aov-error-term-in-r-whats-the-difference-bw-errorid-and-errorid-timevar
# str(litter_N_data)
# hist(litter_N_data$log_N_percent)
litterN.aov <- aov(N_percent~year*treatment+Error(plot/measurement), data = litter_N_data)
summary(litterN.aov)
plot(lm(litterN.aov))
litterN.table<-xtable(litterN.aov)

##########################
# Figure 
##########################
Fig3<-ggplot(data=litter_N_data, aes(y=N_percent, x=treatment, fill=treatment)) +
  stat_boxplot(geom ='errorbar')+  
  geom_boxplot()+
  facet_grid(. ~ year)+
  scale_fill_manual(values = c("Control" = "#FFFFFF", "Low-N" = "#CCCCCC", "High-N" = "#666666"))+
  xlab("Treatment")+
  ylab("Litter N concentration (%)")+
  scale_y_continuous(limit=c(0.2, 0.6))

Fig3<-Fig3+theme_classic() + theme(panel.border = element_blank(), panel.grid.major = element_blank(), 
                                                 axis.line.y = element_line(color="black", size = 0.5, lineend="square"),
                                                 axis.line.x = element_line(color="black", size = 0.5, lineend="square"),
                                                 panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), #sets colors of axes
                                                 plot.title = element_text(hjust=0.05, vjust=-1.8, face="bold", size=22),        #Sets title size, style, location
                                                 axis.title.x=element_text(colour="black", size = 20, vjust=-2),            #sets x axis title size, style, distance from axis #add , face = "bold" if you want bold
                                                 axis.title.y=element_text(colour="black", size = 20, vjust=2),            #sets y axis title size, style, distance from axis #add , face = "bold" if you want bold
                                                 legend.position = c(.9, 0.87),
                                                 # legend.title = element_blank(),
                                                 strip.text = element_text(colour="black", size = 16,face="bold"),                    
                                                 axis.text.x=element_text(colour="black", size = 16, angle=320,hjust=-0.1),                              #sets size and style of labels on axes
                                                 axis.text.y=element_text(colour="black", size = 16),                              #sets size and style of labels on axes
                                                 plot.margin = unit(c(1,3,2,1), "cm"),
                                                 panel.spacing = unit(2, "lines"))
Fig3








############################################################
############################################################
# We also used nested ANOVA to test for effects of N addition on soil ammonium, soil respiration, 
# microbial biomass carbon, qCO2, β-glucosidase, and urease activity. 
# Separate analyses were performed for each parameter, with the four cores nested within the plot 
# in which they were collected. Values for qCO2, soil respiration, MBC, and NH4+ were 
# log-transformed prior to the analyses. 

############################################################
# LOAD AND PROCESS DATA ON soil and microbial properties
############################################################
soils_data<-read.csv("./Data/Resumo_solo.csv", dec=".", header = TRUE, sep = ",", check.names=FALSE )
soils_data$treatment<-recode(soils_data$treatment, 'alto' = "High-N",'baixo' = "Low-N",'controle' = "control")

soils_data$log_beta_g<-log(soils_data$beta_g)
soils_data$log_micro_bmass_C<-log(soils_data$micro_bmass_C)
soils_data$log_nh4<-log(soils_data$NH4)
soils_data$log_urease<-log(soils_data$urease)
soils_data$log_soil_resp<-log(soils_data$soil_resp)
soils_data$log_qco2<-log(soils_data$qco2)

##############################
# ANALYSES
##############################

library(dplyr,tidyr)
# soils_data <- soils_data %>% arrange(treatment)
# soils_data$plot_replicate<-rep(rep(1:5, each=4),length.out=60)
colnames(soils_data)
head(soils_data)
soils_data$plot_replicate<-NA

soils_data = soils_data %>% mutate(plot_replicate = factor(ifelse(plot == "1"|plot == "5"|plot == "2", "1", plot_replicate)))
soils_data = soils_data %>% mutate(plot_replicate = factor(ifelse(plot == "3"|plot == "6"|plot == "8", "2", plot_replicate)))
soils_data = soils_data %>% mutate(plot_replicate = factor(ifelse(plot == "4"|plot == "12"|plot == "9", "3", plot_replicate)))
soils_data = soils_data %>% mutate(plot_replicate = factor(ifelse(plot == "7"|plot == "13"|plot == "11", "4", plot_replicate)))
soils_data = soils_data %>% mutate(plot_replicate = factor(ifelse(plot == "10"|plot == "14"|plot == "15", "5", plot_replicate)))

soils_data$plot_replicate<-as.factor(soils_data$plot_replicate)
soils_data$replicate<-as.factor(soils_data$replicate)

soils_data = soils_data %>% mutate(sample = factor(ifelse(replicate == "a", "1", replicate)))
soils_data = soils_data %>% mutate(sample = factor(ifelse(replicate == "b", "2", replicate)))
soils_data = soils_data %>% mutate(sample = factor(ifelse(replicate == "c", "3", replicate)))
soils_data = soils_data %>% mutate(sample = factor(ifelse(replicate == "d", "4", replicate)))


colnames(soils_data)
soils_data<-soils_data %>% select(number,plot,treatment,plot_replicate,sample,beta_g,
                                  micro_bmass_C,NH4,urease,soil_resp,qco2,log_beta_g,
                                  log_micro_bmass_C,log_nh4,log_urease,log_soil_resp,log_qco2) %>% 
  rename("sampleID"="number","plotID"="plot") %>% 
  arrange(treatment)

head(soils_data)
str(soils_data)

soils_data$plotID<-as.factor(soils_data$plotID)
soils_data$plot_replicate<-as.factor(soils_data$plot_replicate)
soils_data$sample<-as.factor(soils_data$sample)

a<-nlevels(soils_data$treatment)
b<-nlevels(soils_data$plot_replicate)
n<-nlevels(soils_data$sample)

DF_AmongGroups<-(a-1)
DF_AmongGroups
DF_AmongReplicatesWithinGroups<-a*(b-1)
DF_AmongReplicatesWithinGroups
DF_Residual<-a*b*(n-1)
DF_Residual
DF_total<-a*b*n-1
DF_total
a
b
n

# MS_groups<-SS_amonggroups/DF_AmongGroups
# MS_AmongRepsWithinGroups<-SS_replicates_groups/DF_AmongReplicatesWithinGroups
# MS_subsamples<-SS_subsamples/DF_Residual

# F_groups<-MS_groups/MS_AmongRepsWithinGroups  
# F_groups_df_1<-DF_AmongGroups
# F_groups_df_2<-DF_AmongReplicatesWithinGroups

# F_AmongReplicatesWithingGroups<-MS_AmongRepsWithinGroups/# MS_subsamples
# F_AmongReplicatesWithingGroups_df_1<-DF_AmongReplicatesWithinGroups
# F_groups_df_2<-DF_Residual

# https://rcompanion.org/rcompanion/d_07.html
# fit = aov(Protein ~ Tech + Error(Rat), data=Data)
# summary(fit)

####### NH4
head(soils_data)
# fit_log_nh4 = aov(log_nh4 ~ treatment + Error(treatment/plot_replicate), data=soils_data)
fit_log_nh4 = aov(log_nh4 ~ treatment + Error(plotID), data=soils_data)
# to use these see http://online.sfsu.edu/efc/classes/biol458/labs/lab8/Lab-8-R-version.pdf
# fit_log_nh4 = aov(log_nh4 ~ treatment/plotID, data=soils_data)
# fit_log_nh4 = summary(lmer(log_nh4 ~ treatment + (1|plotID), data=soils_data))
# fit_log_nh4 = aov(log_nh4 ~ treatment + plot_replicate:treatment, data=soils_data)
summary(fit_log_nh4)

# res2 <- lm(log_nh4 ~ treatment/plotID, data = soils_data)
# anova(res2)
# fit_log_nh4 = lm(log_nh4 ~ treatment + plot_replicate/sample_replicate, data=soils_data)

# The summary of the aov will produce the correct test for TREATMENT  The test for 
# PLOT can be performed by manually calculating the p-value for the F-test using
# the output for Error:PLOT_REPLICATE and Error:Within.
# Using Mean Sq and Df values to get p-value for H = treatment and Error = Plot
# THIS IS THE SAME AS summary(fit)
# pf(q=.5160/0.3284,
#    df1=2,
#    df2=12,
#    lower.tail=FALSE)

# Using Mean Sq and Df values to get p-value for H = Rat and Error = Within
# # SS of residuals of error plotID/ss of residuals within
pf(q=3.941/2.675,
   df1=12,
   df2=45,
   lower.tail=F)

# pf(q=3.941/2.675,
#    df1=12,
#    df2=30,
#    lower.tail=F)

# MS_groups<-SS_amonggroups/DF_AmongGroups
MS_groups<-1.03/DF_AmongGroups
# MS_AmongRepsWithinGroups<-SS_replicates_groups/DF_AmongReplicatesWithinGroups
MS_AmongRepsWithinGroups<-3.94/DF_AmongReplicatesWithinGroups
# MS_subsamples<-SS_subsamples/DF_Residual
MS_subsamples<-2.68/DF_Residual

# F_groups<-MS_groups/MS_AmongRepsWithinGroups  
F_groups<-MS_groups/MS_AmongRepsWithinGroups  
F_groups_df_1<-DF_AmongGroups
F_groups_df_2<-DF_AmongReplicatesWithinGroups

F_AmongReplicatesWithingGroups<-MS_AmongRepsWithinGroups/MS_subsamples
F_AmongReplicatesWithingGroups_df_1<-DF_AmongReplicatesWithinGroups
F_groups_df_2<-DF_Residual

log_nh4_2.table<-xtable(aov(log_nh4 ~ treatment + Error(plotID), data=soils_data))
log_nh4_2.table<-round(log_nh4_2.table, 2)
write.csv(log_nh4_2.table, file="./Output/log_nh4.aov.table2.csv",na = "")
# plot(fit_log_nh4)
# par(op)
############

head(soils_data)
fit_log_soil_resp2 = aov(log_soil_resp ~ treatment + Error(plotID), data=soils_data)
# fit_log_nh4 = lm(log_nh4 ~ treatment + plot_replicate/sample_replicate, data=soils_data)
summary(fit_log_soil_resp2)

# The summary of the aov will produce the correct test for TREATMENT  The test for 
# PLOT can be performed by manually calculating the p-value for the F-test using
# the output for Error:PLOT_REPLICATE and Error:Within.
# Using Mean Sq and Df values to get p-value for H = treatment and Error = Plot
# THIS IS THE SAME AS summary(fit)
pf(q=0.04/0.11,
   df1=2,
   df2=12,
   lower.tail=FALSE)

# Using Mean Sq and Df values to get p-value for H = Rat and Error = Within
# SS of residuals of error plotID/ss of residuals within
# pf(q=1.3694/3.099,
#    df1=12,
#    df2=30,
#    lower.tail=F)

pf(q=1.3694/3.099,
   df1=12,
   df2=45,
   lower.tail=F)


# MS_groups<-SS_amonggroups/DF_AmongGroups
MS_groups<-0.08/DF_AmongGroups
# MS_AmongRepsWithinGroups<-SS_replicates_groups/DF_AmongReplicatesWithinGroups
MS_AmongRepsWithinGroups<-1.37/DF_AmongReplicatesWithinGroups
# MS_subsamples<-SS_subsamples/DF_Residual
MS_subsamples<-3.1/DF_Residual
MS_subsamples<-3.1/45
# F_groups<-MS_groups/MS_AmongRepsWithinGroups  
F_groups<-MS_groups/MS_AmongRepsWithinGroups  
F_groups_df_1<-DF_AmongGroups
F_groups_df_2<-DF_AmongReplicatesWithinGroups

F_AmongReplicatesWithingGroups<-MS_AmongRepsWithinGroups/MS_subsamples
F_AmongReplicatesWithingGroups_df_1<-DF_AmongReplicatesWithinGroups
F_groups_df_2<-DF_Residual


log_soil_resp2.table<-xtable(aov(log_soil_resp ~ treatment + Error(plotID), data=soils_data))
log_soil_resp2.table<-round(log_soil_resp2.table, 2)
write.csv(log_soil_resp2.table, file="./Output/log_soil_resp.aov.table2.csv",na = "")

############
head(soils_data)
fit_log_micro_bmass_C = aov(log_micro_bmass_C ~ treatment + Error(plotID), data=soils_data)
# fit_log_micro_bmass_C = aov(log_micro_bmass_C ~ treatment + treatment/plotID, data=soils_data)
# fit_log_nh4 = lm(log_nh4 ~ treatment + plot_replicate/sample_replicate, data=soils_data)
summary(fit_log_micro_bmass_C)

# The summary of the aov will produce the correct test for TREATMENT  The test for 
# PLOT can be performed by manually calculating the p-value for the F-test using
# the output for Error:PLOT_REPLICATE and Error:Within.
# Using Mean Sq and Df values to get p-value for H = treatment and Error = Plot
# THIS IS THE SAME AS summary(fit)
# pf(q=.5160/0.3284,
#    df1=2,
#    df2=12,
#    lower.tail=FALSE)

# Using Mean Sq and Df values to get p-value for H = Rat and Error = Within
# SS of residuals of error plotID/ss of residuals within
pf(q=9.805/29.04,
   df1=12,
   df2=45,
   lower.tail=F)
pf(q=9.805/29.04,
   df1=12,
   df2=30,
   lower.tail=F)


# MS_groups<-SS_amonggroups/DF_AmongGroups
MS_groups<-0.36/DF_AmongGroups
# MS_AmongRepsWithinGroups<-SS_replicates_groups/DF_AmongReplicatesWithinGroups
MS_AmongRepsWithinGroups<-9.81/DF_AmongReplicatesWithinGroups
# MS_subsamples<-SS_subsamples/DF_Residual
MS_subsamples<-29.04/DF_Residual
MS_subsamples<-29.04/45
# F_groups<-MS_groups/MS_AmongRepsWithinGroups  
F_groups<-MS_groups/MS_AmongRepsWithinGroups  
F_groups_df_1<-DF_AmongGroups
F_groups_df_2<-DF_AmongReplicatesWithinGroups

F_AmongReplicatesWithingGroups<-MS_AmongRepsWithinGroups/MS_subsamples
F_AmongReplicatesWithingGroups_df_1<-DF_AmongReplicatesWithinGroups
F_groups_df_2<-DF_Residual

log_micro_bmass_C2.table<-xtable(aov(log_micro_bmass_C ~ treatment + Error(plotID), data=soils_data))
log_micro_bmass_C2.table<-round(log_micro_bmass_C2.table, 2)
write.csv(log_micro_bmass_C2.table, file="./Output/log_micro_bmass_C2.aov.table2.csv",na = "")

############
head(soils_data)
fit_log_qco2 = aov(log_qco2 ~ treatment + Error(plotID), data=soils_data)
# + treatment/plotID
# fit_log_nh4 = lm(log_nh4 ~ treatment + plot_replicate/sample_replicate, data=soils_data)
summary(fit_log_qco2)

# The summary of the aov will produce the correct test for TREATMENT  The test for 
# PLOT can be performed by manually calculating the p-value for the F-test using
# the output for Error:PLOT_REPLICATE and Error:Within.
# Using Mean Sq and Df values to get p-value for H = treatment and Error = Plot
# THIS IS THE SAME AS summary(fit)
# pf(q=.5160/0.3284,
#    df1=2,
#    df2=12,
#    lower.tail=FALSE)

# Using Mean Sq and Df values to get p-value for H = Rat and Error = Within
# SS of residuals of error plotID/ss of residuals within
pf(q=12.852/34,
   df1=12,
   df2=45,
   lower.tail=F)
pf(q=12.852/34,
   df1=12,
   df2=30,
   lower.tail=F)


# MS_groups<-SS_amonggroups/DF_AmongGroups
MS_groups<-0.21/DF_AmongGroups
# MS_AmongRepsWithinGroups<-SS_replicates_groups/DF_AmongReplicatesWithinGroups
MS_AmongRepsWithinGroups<-12.85/DF_AmongReplicatesWithinGroups
# MS_subsamples<-SS_subsamples/DF_Residual
MS_subsamples<-34/DF_Residual
MS_subsamples<-34/45
# F_groups<-MS_groups/MS_AmongRepsWithinGroups  
F_groups<-MS_groups/MS_AmongRepsWithinGroups  
F_groups_df_1<-DF_AmongGroups
F_groups_df_2<-DF_AmongReplicatesWithinGroups

F_AmongReplicatesWithingGroups<-MS_AmongRepsWithinGroups/MS_subsamples
F_AmongReplicatesWithingGroups_df_1<-DF_AmongReplicatesWithinGroups
F_groups_df_2<-DF_Residual


log_qco2.table<-xtable(aov(log_qco2 ~ treatment + Error(plotID), data=soils_data))
log_qco2.table<-round(log_qco2.table, 2)
write.csv(log_qco2.table, file="./Output/log_qco22.aov.table2.csv",na = "")

############
head(soils_data)
fit_beta_g = aov(beta_g ~ treatment + Error(plotID), data=soils_data)
# + treatment/plotID
# fit_log_nh4 = lm(log_nh4 ~ treatment + plot_replicate/sample_replicate, data=soils_data)
summary(fit_beta_g)

# The summary of the aov will produce the correct test for TREATMENT  The test for 
# PLOT can be performed by manually calculating the p-value for the F-test using
# the output for Error:PLOT_REPLICATE and Error:Within.
# Using Mean Sq and Df values to get p-value for H = treatment and Error = Plot
# THIS IS THE SAME AS summary(fit)
# pf(q=.5160/0.3284,
#    df1=2,
#    df2=12,
#    lower.tail=FALSE)

# Using Mean Sq and Df values to get p-value for H = Rat and Error = Within
# SS of residuals of error plotID/ss of residuals within
pf(q=8118/32218,
   df1=12,
   df2=45,
   lower.tail=F)



# MS_groups<-SS_amonggroups/DF_AmongGroups
MS_groups<-1646/DF_AmongGroups
# MS_AmongRepsWithinGroups<-SS_replicates_groups/DF_AmongReplicatesWithinGroups
MS_AmongRepsWithinGroups<-8118/DF_AmongReplicatesWithinGroups
# MS_subsamples<-SS_subsamples/DF_Residual
MS_subsamples<-32218/DF_Residual
MS_subsamples<-32218/45
# F_groups<-MS_groups/MS_AmongRepsWithinGroups  
F_groups<-MS_groups/MS_AmongRepsWithinGroups  
F_groups_df_1<-DF_AmongGroups
F_groups_df_2<-DF_AmongReplicatesWithinGroups

F_AmongReplicatesWithingGroups<-MS_AmongRepsWithinGroups/MS_subsamples
F_AmongReplicatesWithingGroups_df_1<-DF_AmongReplicatesWithinGroups
F_groups_df_2<-DF_Residual


beta_g.table<-xtable(aov(beta_g ~ treatment + Error(plotID), data=soils_data))
beta_g.table<-round(beta_g.table, 2)
write.csv(beta_g.table, file="./Output/beta_g2.aov.table2.csv",na = "")

############
head(soils_data)
fit_urease = aov(urease ~ treatment + Error(plotID), data=soils_data)
# + treatment/plotID
# fit_log_nh4 = lm(log_nh4 ~ treatment + plot_replicate/sample_replicate, data=soils_data)
summary(fit_urease)
# The summary of the aov will produce the correct test for TREATMENT  The test for 
# PLOT can be performed by manually calculating the p-value for the F-test using
# the output for Error:PLOT_REPLICATE and Error:Within.
# Using Mean Sq and Df values to get p-value for H = treatment and Error = Plot
# THIS IS THE SAME AS summary(fit)
pf(q=24471/253973,
   df1=2,
   df2=12,
   lower.tail=FALSE)

# Using Mean Sq and Df values to get p-value for H = Rat and Error = Within
# SS of residuals of error plotID/ss of residuals within
pf(q=253973/22057,
   df1=12,
   df2=45,
   lower.tail=F)

# pf(q=253973/22057,
#    df1=12,
#    df2=30,
#    lower.tail=F)

# MS_groups<-SS_amonggroups/DF_AmongGroups
MS_groups<-24471/DF_AmongGroups
# MS_AmongRepsWithinGroups<-SS_replicates_groups/DF_AmongReplicatesWithinGroups
MS_AmongRepsWithinGroups<-253973/DF_AmongReplicatesWithinGroups
# MS_subsamples<-SS_subsamples/DF_Residual
MS_subsamples<-22057/DF_Residual
MS_subsamples<-22057/45
# F_groups<-MS_groups/MS_AmongRepsWithinGroups  
F_groups<-MS_groups/MS_AmongRepsWithinGroups  
F_groups_df_1<-DF_AmongGroups
F_groups_df_2<-DF_AmongReplicatesWithinGroups

F_AmongReplicatesWithingGroups<-MS_AmongRepsWithinGroups/MS_subsamples
F_AmongReplicatesWithingGroups_df_1<-DF_AmongReplicatesWithinGroups
F_groups_df_2<-DF_Residual

urease.table<-xtable(aov(urease ~ treatment + Error(plotID), data=soils_data))
urease.table<-round(urease.table, 2)
write.csv(urease.table, file="./Output/urease2.aov.table2.csv",na = "")

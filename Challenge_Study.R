#This is the code that was used for challenge study analysis as a part of the
#manuscript ": Listeria monocytogenes, Salmonella enterica, 
#and enterohemorrhagic Escherichia coli show significant die-off in both conventional 
#and organic controlled environment agriculture nutrient solution over a five-day period."

## install if needed (do this exactly once):
## install.packages("usethis")

library(usethis)
use_git_config(user.name = "tamardigrade", user.email = "tmwalsky@gmail.com")

#Load  packages
library("tidyverse")
library(lmerTest)
library(emmeans)
library(nlsMicrobio)
library(grid)
library(gridExtra)
library(extrafont)
library(broom.mixed)
font_import() # Make sure to type 'Y' when prompted
loadfonts(device = "win", quiet = TRUE)

#save the R environment
#save.image(file='C:/Users/tmwal/OneDrive/Desktop/Challengestudydata04242026.RData')

#See the following link for help with interpreting linear models, regression and anova
#https://cscu.cornell.edu/workshop/interpreting-linear-models-regression-and-anova/

#load saved R environment
#load("C:/Users/tmwal/OneDrive/Desktop/Challengestudydata04242026.RData")
load("/local1/workdir1/tw488/CEA_RData_RScript/Challengestudydata04242026.RData")
save.image(file='/local1/workdir1/tw488/CEA_RData_RScript/Challengestudydata04242026.RData')


#import Total Plate Count (TPC) files 
TPC<-read.csv('/local1/workdir1/tw488/CEA_RData_RScript/Challenge_study/Induvidual_Water_Samples_for_manuscipt_07022025.csv')
TPC<-TPC%>%na.omit()
TPC$Pond_Type_Number <- paste(TPC$Pond_Type,TPC$Pond_Number)

#calculate TPC average values (using replicates)
TPC1 <- TPC %>% group_by(Sample_Name) %>% mutate(Average_Conc=mean(Concentration))
TPC1 <- TPC1 %>% group_by(Sample_Name) %>% mutate(Log_Average_Conc=log10(Average_Conc))%>%ungroup

#Use lmer to make mixed model to test impact of pond type on log count data using unaveraged replicates
  # Pond Number:fixed effect
  # Location:fixed effect
TPC2<-TPC1%>%distinct(Pond_Type_Number,Pond_Type,Pond_Number, Log_Average_Conc, .keep_all = TRUE)
TPC2$Pond_Type_Number <- as.factor(TPC2$Pond_Type_Number)


#TPC_lmer <- lme4::lmer(log10(Average_Conc) ~ Pond_Type+ Location+(1|Pond_Type_Number), TPC2)
TPC_lmer <- lm(log10(Average_Conc) ~ Pond_Type+ Location, TPC2)

summary(TPC_lmer)

TPC_lmer_df <- tidy(TPC_lmer, conf.int = TRUE)
write.csv(TPC_lmer_df, "/local1/workdir1/tw488/CEA_RData_RScript/Challenge_study/APC_lm_emmeans/TPC_lmer_df.csv", row.names = FALSE)


#Figure out effect size for pond type
emmeans_TPC_lmer<-emmeans(TPC_lmer, pairwise ~ Pond_Type, type = "response")
emmeans_TPC_lmer

#https://cran.r-project.org/web/packages/emmeans/vignettes/transformations.html
#Alternative scale section for percent differences instead of ratios

pct.diff.tran10 <- list(
  linkfun = function(mu) log10(mu/100 + 1),
  #these next 3 lines are changed to reflect log10, not natural log
  linkinv = function(eta) 100 * (10^(eta) - 1),
  mu.eta  = function(eta) 100 * log(10) * 10^(eta),
  name    = "log10(pct.diff)"
)

emmeans_TPC_lmer_df <- as.data.frame(update(pairs(emmeans_TPC_lmer, type = "response"), 
       tran = pct.diff.tran10, inv.lbl = "pct.diff", adjust = "none",
       infer = c(TRUE, TRUE)))
emmeans_TPC_lmer_df

#Figure out effect size for location
emmeans_TPC_lmer_Loc<-emmeans(TPC_lmer, pairwise ~ Location, type = "response")
emmeans_TPC_lmer_Loc


pct.diff.tran10 <- list(
  linkfun = function(mu) log10(mu/100 + 1),
  #these next 3 lines are changed to reflect log10, not natural log
  linkinv = function(eta) 100 * (10^(eta) - 1),
  mu.eta  = function(eta) 100 * log(10) * 10^(eta),
  name    = "log10(pct.diff)"
)

emmeans_TPC_lmer_Loc_df <- as.data.frame(update(pairs(emmeans_TPC_lmer_Loc, type = "response"), 
       tran = pct.diff.tran10, inv.lbl = "pct.diff", adjust = "none",
       infer = c(TRUE, TRUE)))
emmeans_TPC_lmer_Loc_df

#Save APC emmeans details to 
write.csv(emmeans_TPC_lmer_Loc_df, "/local1/workdir1/tw488/CEA_RData_RScript/Challenge_study/APC_lm_emmeans/emmeans_TPC_lmer_Loc_df.csv", row.names = FALSE)
write.csv(emmeans_TPC_lmer_df, "/local1/workdir1/tw488/CEA_RData_RScript/Challenge_study/APC_lm_emmeans/emmeans_TPC_lmer_df.csv", row.names = FALSE)


# Simple Linear Model
TPC_final <- lm(Log_Average_Conc ~ Pond_Type + Location, data = TPC2)

# Check diagnostic plots 
plot(TPC_final)
anova(TPC_final)
summary(TPC_final)


#Calculate the average log concentration of conventional vs organic samples
#TPC3 <- TPC2 %>% group_by(Pond_Type) %>% 
 # mutate(Average_Conc_byPondType=mean(Average_Conc))%>%
  #mutate(Log_Average_Conc_byPondType=log10(Average_Conc_byPondType))%>%
  #mutate(log_SD_Average_Conc_byPond=sd(Log_Average_Conc)) %>%ungroup
TPC3v2 <- TPC2 %>% group_by(Pond_Type) %>% 
 # mutate(Log_Conc=log10(Average_Conc))%>%
  mutate(Log_Average_Conc_byPondType=mean(Log_Average_Conc))%>%
  mutate(log_SD_Average_Conc_byPond=sd(Log_Average_Conc)) %>%ungroup

#TPC_logavgconcbypondtype<-TPC3%>%distinct(Pond_Type, Log_Average_Conc_byPondType, .keep_all = TRUE)
TPC_logavgconcbypondtypev2<-TPC3v2%>%distinct(Pond_Type, Log_Average_Conc_byPondType,log_SD_Average_Conc_byPond, .keep_all = TRUE)

write.csv(TPC_logavgconcbypondtypev2, "/local1/workdir1/tw488/CEA_RData_RScript/Challenge_study/PLate_Count_Values_for_Paper/TPC_logavgconcbypondtypev2.csv", row.names = FALSE)



#INOCULATION DATA
#import pathogen data
Day0<-read.csv("/local1/workdir1/tw488/CEA_RData_RScript/Challenge_study/DAY0_COMBINED_FOR_R_updated_v2.csv")
Day1<-read.csv("/local1/workdir1/tw488/CEA_RData_RScript/Challenge_study/DAY1_COMBINED_FOR_R_updated_v2.csv")
Day2<-read.csv("/local1/workdir1/tw488/CEA_RData_RScript/Challenge_study/DAY2_COMBINED_FOR_R_updated_v2.csv")
Day5<-read.csv("/local1/workdir1/tw488/CEA_RData_RScript/Challenge_study/DAY5_combined_for_R_updated_v2.csv")

#MERGE data all in one frame
df <- rbind(Day0,Day1) 
df <- rbind(df,Day2) 
df <- rbind(df,Day5) 

#add new groupings/ sample name combos
df$Sample_ID_Aliquot<-paste(df$Day,df$Pond_Type_Number,df$Aliquot,sep = '_')
df$Sample_ID<-paste(df$Day,df$Pond_Type_Number,sep = '_')
df$Pond_Type_Number<-as.character(df$Pond_Type_Number)

#Listeria monocytogenes______________________________________________________

#sort out Lmono
Lmono <- df %>% filter(Bacteria=="LMONO")

#Sort out Lmono NC, check that it looks ok
LmonoNC <- df %>% filter(Bacteria=="LMONO")%>% filter(Aliquot=="NC")

#Take out Lmono NC, calc average values
Lmono_A_C <- df %>% filter(Bacteria=="LMONO")%>% filter(Aliquot!="NC")


#find limits of detection (LOD) for undiluted and 0.1 dilution plates
#detection limit (CFU/mL) = (# colonies ÷(the amount plated on the agar plate)) × dilution factor 
LOD_CFU_mL_undiluted = (1/ 0.05)*1 
LOD_CFU_mL_0.1 = (1/ 0.05)*10
log10(LOD_CFU_mL_undiluted)
log10(LOD_CFU_mL_0.1)

#Calculate LOD times * 0.25 for imputation
Impute_LOD_CFU_mL_0.1=LOD_CFU_mL_0.1*0.25
Impute_LOD_CFU_mL_undiluted=LOD_CFU_mL_undiluted*0.25
log10(Impute_LOD_CFU_mL_0.1)
log10(Impute_LOD_CFU_mL_undiluted)


#Replace 0 in Avg_Conc column of Day 5 with LOD_CFU_mL_undiluted (the Limit of detection for undiluted plates)
#Replace 0 in Avg_Conc column of Day 5 with LOD_CFU_mL_0.1 (the Limit of detection for 0.1 dilution plates)

Lmono_A_C <- Lmono_A_C %>%
  mutate(
    Final_Concentration = as.numeric(Final_Concentration),
    Below_LOD = Final_Concentration == 0,          
    LOD = case_when(
      Day == "DAY2" ~ Impute_LOD_CFU_mL_0.1,
      Day == "DAY5" ~ Impute_LOD_CFU_mL_undiluted,
      TRUE          ~ NA_real_                     #no nondetects on DAY0/DAY1
    ),
    Final_Concentration = ifelse(Below_LOD, LOD, Final_Concentration)
  )

#Check: this should return nothing
Lmono_A_C %>% filter(Below_LOD & is.na(LOD))

#calculate average values by replicate
Lmono_A_C_byREP <- Lmono_A_C %>% 
  group_by(Day,Pond_Type_Number, Aliquot, Dilution) %>% 
  mutate(Average_Conc=mean(Final_Concentration)) %>% 
  mutate(Log_Average_Conc=log10(Average_Conc))%>%
  ungroup

Lmono_A_C_1 <- Lmono_A_C_byREP %>%distinct(Day,Pond_Type_Number, Aliquot, .keep_all = TRUE)


#plot

mycolslmono <- rep(c("#E69F00", 
                    "#D55E00", 
                    "#56B4E9",
                    "#0072B2"), 
                  length.out = length(unique(Lmono_A_C_1$Sample_ID)))

Lmono_A_C_1$Day <- factor(Lmono_A_C_1$Day, levels = c("DAY0", "DAY1", "DAY2", "DAY5"), 
                                 labels = c("Day 0", "Day 1", "Day 2","Day 5"))

y.expression <- expression(Concentration~ (log [10] ~CFU/mL))

#Make a continuous day column

Lmono_A_C_1<-Lmono_A_C_1%>%mutate(Day_numeric=Day)%>%mutate(Day_numeric = str_remove_all(Day, "[^0-9.-]")) %>%
  mutate(Day_numeric = as.numeric(Day_numeric))

#plot
r5 <- ggplot(Lmono_A_C_1, aes(x = Day_numeric, y = Log_Average_Conc, 
                                         shape = factor(Pond_Type_Number), 
                                         color = factor(Pond_Type_Number))) + 
  # Map both color and shape to the same factor to merge the legend
  geom_point(size = 12, alpha = 0.9) + 
  labs(title = expression( ~ italic("Listeria monocytogenes")),
       x = "Day", 
       y = y.expression) +  
  theme_bw() +
  # Use show.legend = FALSE to keep the dashed line out of your point legend
  geom_hline(yintercept = 1.3, linetype = "dashed", show.legend = FALSE) +
  ylim(1, 5) +
  # Ensure the manual color scale name matches the shape scale (both are NULL/blank)
  scale_color_manual(values = mycolslmono) +
  # Optional: force shape values if you don't like the defaults
  scale_shape_manual(values = c(16, 17, 15, 18)) + 
  theme(text = element_text(family = "Times New Roman"),
        strip.text = element_text(size = 40, face = "bold"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        legend.title = element_blank(), # This helps merging
        legend.text = element_text(size = 30),
        axis.text.x = element_text(angle = 0, vjust = 0.3, hjust = 1, size = 40, margin = margin(0,0,20,0)),
        axis.text.y = element_text(angle = 0, vjust = 0.3, hjust = 0.5, size = 40, margin = margin(0,0,0,20)),
        axis.title.y = element_text(size = 50),
        axis.title.x = element_text(size = 50),
        plot.title = element_text(size = 55, face = "italic", hjust = 0.5))

r5

ggsave(
  filename = "/local1/workdir1/tw488/CEA_RData_RScript/Challenge_study/Lmonolin_contxaxis_09292026.tiff", 
  plot = last_plot(),       
  device = "tiff",          
  units = "in", 
  width = 25, 
  height = 15, 
  dpi = 700, 
  compression = "lzw",
  type = "cairo"            
)

#COMPARING LINEAR AND MAFART WEIBULL MODELS
# Fit a linear model
Lmonomodel <- lm(Log_Average_Conc ~ Day_numeric, data = Lmono_A_C_1)

# Generate 4 diagnostic plots (Residuals vs Fitted is first)
plot(Lmonomodel) 
BIC(Lmonomodel)
summary(Lmonomodel)

# Get 95% Confidence Intervals
confint(Lmonomodel)
Lmonomodelparameters_Lm <- summary(Lmonomodel)$coefficients
print(Lmonomodelparameters_Lm)

#Extract the slope (coefficient of the 'time' variable)
slope <- coef(Lmonomodel)["Day_numeric"]

# Calculate the D-value (-1 / slope)
d_value = -1 / slope
d_value

# Calculate the confidence intervals
ci_slope <- confint(Lmonomodel)["Day_numeric", ]
d_ci <- rev(-1 / ci_slope)    
d_ci

# To see only the Residuals vs Fitted plot:
plot(Lmonomodel, which = 1)
ggsave("/local1/workdir1/tw488/CEA_RData_RScript/Challenge_study/Check_linearity/Lmonoresiduals.tiff", units="in", width=20, height=15, dpi=500, compression = 'lzw')

# qqplotplot:
plot(Lmonomodel, which = 2)
ggsave("/local1/workdir1/tw488/CEA_RData_RScript/Challenge_study/Check_linearity/Lmonoqq.tiff", units="in", width=20, height=15, dpi=500, compression = 'lzw')

# Fit the Lmono data to a weibull distribution (mafart)
mafart_weibull_Lm <- nls(Log_Average_Conc ~ log10_N0 - (Day_numeric/delta)^p, 
                        data =Lmono_A_C_1, 
                        start = list(log10_N0 = 4.661987, delta = 1, p = 0.5))


overview(mafart_weibull_Lm)
plotfit(mafart_weibull_Lm, smooth = TRUE)
BIC(mafart_weibull_Lm)
mafartparameters_Lm <- summary(mafart_weibull_Lm)$coefficients
print(mafartparameters_Lm)

#check residuals
plot(fitted(mafart_weibull_Lm), residuals(mafart_weibull_Lm))
abline(h = 0, lty = 2, col = "red")

#Use lmer to make mixed model to test impact of pond type on log count data using unaveraged aliquot values
# use log(Average_Conc) in the formula to inform future emmeans code that this data is log transformed
# see this vignette for a full explanation: https://cran.r-project.org/web/packages/emmeans/vignettes/transformations.html

class(Lmono_A_C_1)
Lmono_lmer <- lmerTest::lmer((Log_Average_Conc) ~ Day * Pond_Type+(1|Pond_Type_Number), Lmono_A_C_1)
summary(Lmono_lmer)
anova(Lmono_lmer)

#add ',type="response" to this emmeans code so it recognizes
#the values are log10 transformed and the responses emmeans provides are exponentiated
#so you dont have to exponentiate by hand
Lmono_lmer_emmeans <- emmeans(Lmono_lmer, consec ~ Day,  type = "response", lmer.df = "satterthwaite")


#To calculate % change 
Lmono_lmer_emmeans_df <- as.data.frame(update(
  contrast(Lmono_lmer_emmeans, method = "consec"),   # DAY1/DAY0, DAY2/DAY1, DAY5/DAY2
  tran = pct.diff.tran10, inv.lbl = "pct.diff",
  adjust = "none",                                   # or remove to use the mvt default
  infer = c(TRUE, TRUE)))

Lmono_lmer_emmeans_df

#Save Lmono_lmer_emmeans details to 
write.csv(Lmono_lmer_emmeans_df, "/local1/workdir1/tw488/CEA_RData_RScript/Challenge_study/Pathogen_lm_emmeans/Lmono_lmer_emmeans.csv", row.names = FALSE)


#Average the aliquots together to find overall average values by organic/ conventional
Lmono_A_C_bypond <- Lmono_A_C_1%>%
  group_by(Day,Pond_Type) %>% 
  mutate(Average_Log_Conc_byPond=mean(Log_Average_Conc))%>%
  mutate(log_SD_Average_Conc_byPond=sd(Log_Average_Conc)) %>%
  mutate(log_SER_Average_Conc_byPond=(sd(Log_Average_Conc))/(sqrt(3))) %>%
  ungroup

Lmono_A_C_bypond <- Lmono_A_C_bypond %>%distinct(Day,Pond_Type, .keep_all = TRUE)

Lmono_A_C_bypond_fortable <-Lmono_A_C_bypond%>%select(Sample_ID, Bacteria,Average_Log_Conc_byPond,log_SD_Average_Conc_byPond)%>%
  mutate(Average_Log_Conc_byPond = round(Average_Log_Conc_byPond, digits = 2))%>%
  mutate(log_SD_Average_Conc_byPond = round(log_SD_Average_Conc_byPond, digits = 2))
  
#SALMONELLA________________________________________________________________

#sort out Salm
Salm <- df %>% filter(Bacteria=="SALM")
#sort out Salm NC, check that it looks ok
SalmNC <- df %>% filter(Bacteria=="SALM")%>% filter(Aliquot=="NC")
#take out Salm NC, calc avg values
Salm_A_C <- df %>% filter(Bacteria=="SALM")%>% filter(Aliquot!="NC")

#calculate average values by replicate
Salm_A_C_byREP <- Salm_A_C %>% 
  group_by(Day,Pond_Type_Number, Aliquot, Dilution) %>% 
  mutate(Average_Conc=mean(Final_Concentration)) %>% 
  mutate(Log_Average_Conc=log10(Average_Conc))%>%
  ungroup

Salm_A_C_1 <- Salm_A_C_byREP %>%distinct(Day,Pond_Type_Number, Aliquot, .keep_all = TRUE)
str(Salm_A_C_1)
Salm_A_C_1$Day<-as.factor(Salm_A_C_1$Day)


#Make a continuous day column

Salm_A_C_1<-Salm_A_C_1%>%mutate(Day_numeric=Day)%>%mutate(Day_numeric = str_remove_all(Day, "[^0-9.-]")) %>%
  mutate(Day_numeric = as.numeric(Day_numeric))

mycolsSalm <- rep(c("#E69F00", 
                    "#D55E00", 
                    "#56B4E9",
                    "#0072B2"), 
                  length.out = length(unique(Salm_A_C_1$Sample_ID)))

Salm_A_C_1$Day <- factor(Salm_A_C_1$Day, levels = c("DAY0", "DAY1", "DAY2", "DAY5"), 
                         labels = c("Day 0", "Day 1", "Day 2","Day 5"))

rsalm <- ggplot(Salm_A_C_1, aes(x = Day_numeric, y = Log_Average_Conc, 
                                         shape = factor(Pond_Type_Number), 
                                         color = factor(Pond_Type_Number))) + 
  # Map both color and shape to the same factor to merge the legend
  geom_point(size = 12, alpha = 0.9) + 
  labs(title = expression(~italic("Salmonella enterica")),
       x = "Day", 
       y = y.expression) +  
  theme_bw() +
  # Use show.legend = FALSE to keep the dashed line out of your point legend
  geom_hline(yintercept = 1.3, linetype = "dashed", show.legend = FALSE) +
  ylim(1, 5) +
  # Ensure the manual color scale name matches the shape scale (both are NULL/blank)
  scale_color_manual(values = mycolsSalm) +
  # Optional: force shape values if you don't like the defaults
  scale_shape_manual(values = c(16, 17, 15, 18)) + 
  theme(text = element_text(family = "Times New Roman"),
        strip.text = element_text(size = 40, face = "bold"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        legend.title = element_blank(),
        legend.text = element_text(size = 30),
        axis.text.x = element_text(angle = 0, vjust = 0.3, hjust = 1, size = 40, margin = margin(0,0,20,0)),
        axis.text.y = element_text(angle = 0, vjust = 0.3, hjust = 0.5, size = 40, margin = margin(0,0,0,20)),
        axis.title.y = element_text(size = 50),
        axis.title.x = element_text(size = 50),
        plot.title = element_text(size = 55, face = "italic", hjust = 0.5))

rsalm

ggsave(
  filename = "/local1/workdir1/tw488/CEA_RData_RScript/Challenge_study/Salm_10122026.tiff", 
  plot = last_plot(),     
  device = "tiff",          
  units = "in", 
  width = 25, 
  height = 15, 
  dpi = 700, 
  compression = "lzw",
  type = "cairo")

#COMPARING LINEAR AND MAFART WEIBULL MODELS

# Fit a model
Salmmodel <- lm(Log_Average_Conc ~ Day_numeric, data = Salm_A_C_1)
BIC(Salmmodel)

# To see only the Residuals vs Fitted plot:
plot(Salmmodel, which = 1)
ggsave("/local1/workdir1/tw488/CEA_RData_RScript/Challenge_study/Check_linearity/Salmresiduals.tiff", units="in", width=20, height=15, dpi=500, compression = 'lzw')

# qqplotplot:
plot(Salmmodel, which = 2)
ggsave("/local1/workdir1/tw488/CEA_RData_RScript/Challenge_study/Check_linearity/Salmqq.tiff", units="in", width=20, height=15, dpi=500, compression = 'lzw')


#divide Salm_A_C_1 in to organic and conventional
Salm_A_C_1_O<-Salm_A_C_1%>%filter(Pond_Type=="Organic")
Salm_A_C_1_C<-Salm_A_C_1%>%filter(Pond_Type=="Conventional")

#linear model
Salmmodel_C <- lm(Log_Average_Conc ~ Day_numeric, data = Salm_A_C_1_C)
BIC(Salmmodel_C)
Salmmodel_O <- lm(Log_Average_Conc ~ Day_numeric, data = Salm_A_C_1_O)
BIC(Salmmodel_O)

#mafart conventional
mafart_weibull_SC <- nls(Log_Average_Conc ~ log10_N0 - (Day_numeric/delta)^p, 
                      data =Salm_A_C_1_C, 
                      start = list(log10_N0 = 4.556303, delta = 2, p = 1))


overview(mafart_weibull_SC)
plotfit(mafart_weibull_SC, smooth = TRUE)
BIC(mafart_weibull_SC)
mafartparameters_SC <- summary(mafart_weibull_SC)$coefficients
print(mafartparameters_SC)

plot(fitted(mafart_weibull_SC), residuals(mafart_weibull_SC))
abline(h = 0, lty = 2, col = "red")

#mafart ORGANIC
mafart_weibull_SO <- nls(Log_Average_Conc ~ log10_N0 - (Day_numeric/delta)^p, 
                         data =Salm_A_C_1_O, 
                         start = list(log10_N0 = 4.556303, delta = 2, p = 1))


overview(mafart_weibull_SO)
plotfit(mafart_weibull_SO, smooth = TRUE)
BIC(mafart_weibull_SO)
mafartparameters_SO <- summary(mafart_weibull_SO)$coefficients
print(mafartparameters_SO)

plot(fitted(mafart_weibull_SO), residuals(mafart_weibull_SO))
abline(h = 0, lty = 2, col = "red")

#Use lmer to make mixed model to test impact of pond type on log count data using unaveraged aliquots

Salm_lmer <- lmerTest::lmer(log10(Average_Conc) ~ Day * Pond_Type+(1|Pond_Type_Number), Salm_A_C_1)
summary(Salm_lmer)
anova(Salm_lmer)

#add ',type="response" to this emmeans code so it recognizes
#the values are log10 transformed and the responses emmeans provides are exponentiated
#so you dont have to exponentiate by hand

Salm_lmer_emmeans <- emmeans(Salm_lmer, consec ~ Day| Pond_Type,type="response", lmer.df = "satterthwaite") 


#To calculate % change 

Salm_lmer_emmeans_df <- as.data.frame(update(
  contrast(Salm_lmer_emmeans, method = "consec"),   # DAY1/DAY0, DAY2/DAY1, DAY5/DAY2
  tran = pct.diff.tran10, inv.lbl = "pct.diff",
  adjust = "none",                                   # or remove to use the mvt default
  infer = c(TRUE, TRUE)))

Salm_lmer_emmeans_df

#Save Salm_lmer_emmeans details to 
write.csv(Salm_lmer_emmeans_df, "/local1/workdir1/tw488/CEA_RData_RScript/Challenge_study/Pathogen_lm_emmeans/Salm_lmer_emmeans.csv", row.names = FALSE)

#calculate averages

Salm_A_C_bypond <- Salm_A_C_1%>%
  group_by(Day,Pond_Type) %>% 
  mutate(Average_Log_Conc_byPond=(mean(Log_Average_Conc)))%>%
  mutate(log_SD_Average_Conc_byPond=sd(Log_Average_Conc)) %>%
  mutate(log_SER_Average_Conc_byPond=(sd(Log_Average_Conc))/(sqrt(3))) %>%
  ungroup


Salm_A_C_bypond <- Salm_A_C_bypond %>%distinct(Day,Pond_Type, .keep_all = TRUE)

Salm_A_C_bypond_fortable <-Salm_A_C_bypond%>%select(Sample_ID, Bacteria,Average_Log_Conc_byPond,log_SD_Average_Conc_byPond)%>%
  mutate(Average_Log_Conc_byPond = round(Average_Log_Conc_byPond, digits = 2))%>%
  mutate(log_SD_Average_Conc_byPond = round(log_SD_Average_Conc_byPond, digits = 2))


#ECOLI_________________________________________________________________

#import E.coli data with select dilutions
Ecoli_select<-read.csv('/local1/workdir1/tw488/CEA_RData_RScript/Challenge_study/Ecoli_Data_by_Day_forR_07082025_v2.csv')
Ecoli_select$Pond_Type_Number <- paste(Ecoli_select$Pond_Type,Ecoli_select$Pond_Number)

#Make sample ID column
Ecoli_select$Sample_ID_Aliquot<-paste(Ecoli_select$Day,Ecoli_select$Pond_Type_Number,Ecoli_select$Aliquot,sep = '_')
Ecoli_select$Sample_ID<-paste(Ecoli_select$Day,Ecoli_select$Pond_Type_Number,sep = '_')

#calculate average concentration and log average concentration
Ecoli_select1 <- Ecoli_select %>% group_by(Day,Pond_Type_Number,Aliquot,Dilution,Date) %>% 
  mutate(Average_Conc=mean(Final_Conc))%>% 
  mutate(Log_Average_Conc=log10(Average_Conc))%>% ungroup

#get rid of duplicates 
Ecoli_select1<-Ecoli_select1 %>% distinct(Day,Pond_Type_Number,Aliquot,Dilution,Date, .keep_all = TRUE)

#convert sample to character from integer
Ecoli_select1$Pond_Type<-as.character(Ecoli_select1$Pond_Type)
str(Ecoli_select1)

#remove NC samples from Ecoli_select1 to use in mixed model
Ecoli_select1_noNC<-Ecoli_select1%>%  filter(Aliquot!="NC")

##Use lmer to make mixed model to test impact of pond type on log count data using unaveraged aliquots
Ecoli_lmer <- lmerTest::lmer(log10(Average_Conc) ~ Day * Pond_Type+(1|Pond_Type_Number), Ecoli_select1_noNC)
summary(Ecoli_lmer)
anova(Ecoli_lmer)

#add ',type="response" to this emmeans code so it recognizes
#the values are log10 transformed and the responses emmeans provides are exponentiated
#so you dont have to exponentiate by hand

Ecoli_lmer_emmeans <- emmeans(Ecoli_lmer, consec ~ Day| Pond_Type,  type = "response", lmer.df = "satterthwaite")

#To calculate % change 
Ecoli_lmer_emmeans_df <- as.data.frame(update(
  contrast(Ecoli_lmer_emmeans, method = "consec"),   # DAY1/DAY0, DAY2/DAY1, DAY5/DAY2
  tran = pct.diff.tran10, inv.lbl = "pct.diff",
  adjust = "none",                                   # or remove to use the mvt default
  infer = c(TRUE, TRUE)))

Ecoli_lmer_emmeans_df

#Save Lmono_lmer_emmeans details to 
write.csv(Ecoli_lmer_emmeans_df, "/local1/workdir1/tw488/CEA_RData_RScript/Challenge_study/Pathogen_lm_emmeans/Ecoli_lmer_emmeans.csv", row.names = FALSE)


#remove NC samples, adjust the Sample ID so it has NC attached and then merge with Ecoli_select1_noNC df
Ecoli_select1_NC <- Ecoli_select1%>%
  filter(Aliquot=="NC")
Ecoli_select1_NC$Sample_ID<-paste(Ecoli_select1_NC$Day,Ecoli_select1_NC$Pond_Type_Number,Ecoli_select1_NC$Aliquot,sep = '_')

#Combine data frames
Ecoli_select1v2<-merge(Ecoli_select1_noNC,Ecoli_select1_NC, all=TRUE)

#Arrange so its sorted properly each time
Ecoli_select1v2<-Ecoli_select1v2%>%arrange(Sample_ID_Aliquot)


mycolsEcoli <- rep(c("#E69F00", 
                     "#D55E00", 
                     "#56B4E9",
                     "#0072B2"), 
              length.out = length(unique(Ecoli_select1_noNC$Sample_ID)))

Ecoli_select1_noNC$Day <- factor(Ecoli_select1_noNC$Day, levels = c("DAY0", "DAY1", "DAY2", "DAY5"), 
                         labels = c("Day 0", "Day 1", "Day 2","Day 5"))

Ecoli_select1_noNC<-Ecoli_select1_noNC%>%mutate(Day_numeric=Day)%>%mutate(Day_numeric = str_remove_all(Day, "[^0-9.-]")) %>%
  mutate(Day_numeric = as.numeric(Day_numeric))

LOD.expression <- expression("The Limit of Detection is 1.3"~log [10] ~CFU/mL)

recoli <- ggplot(Ecoli_select1_noNC, aes(x = Day_numeric, y = Log_Average_Conc, 
                                         shape = factor(Pond_Type_Number), 
                                         color = factor(Pond_Type_Number))) + 
  geom_point(size = 12, alpha = 0.9) + 
  labs(title = expression("Enterohemorrhagic" ~ italic("E. coli")),
       x = "Day", 
       y = y.expression) +  
  theme_bw() +
  # Use show.legend = FALSE to keep the dashed line out of your point legend
  geom_hline(yintercept = 1.3, linetype = "dashed", show.legend = FALSE) +
  ylim(1, 5) +
  # Ensure the manual color scale name matches the shape scale (both are NULL/blank)
  scale_color_manual(values = mycolsEcoli) +
  # Optional: force shape values if you don't like the defaults
  scale_shape_manual(values = c(16, 17, 15, 18)) + 
  theme(text = element_text(family = "Times New Roman"),
        strip.text = element_text(size = 40, face = "bold"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        legend.title = element_blank(), # This helps merging
        legend.text = element_text(size = 30),
        axis.text.x = element_text(angle = 0, vjust = 0.3, hjust = 1, size = 40, margin = margin(0,0,20,0)),
        axis.text.y = element_text(angle = 0, vjust = 0.3, hjust = 0.5, size = 40, margin = margin(0,0,0,20)),
        axis.title.y = element_text(size = 50),
        axis.title.x = element_text(size = 50),
        plot.title = element_text(size = 55, face = "italic", hjust = 0.5))

recoli


ggsave(
  filename = "/local1/workdir1/tw488/CEA_RData_RScript/Challenge_study/Ecoli_10122026.tiff", 
  plot = last_plot(),       
  device = "tiff",          
  units = "in", 
  width = 25, 
  height = 15, 
  dpi = 700, 
  compression = "lzw",
  type = "cairo")

#COMPARING LINEAR AND MAFART WEIBULL MODELS
#make a continuous day column

Ecoli_select1_noNC<-Ecoli_select1_noNC%>%mutate(Day_numeric=Day)%>%mutate(Day_numeric = str_remove_all(Day, "[^0-9.-]")) %>%
  mutate(Day_numeric = as.numeric(Day_numeric))

# Fit a model
Ecolimodel <- lm(Log_Average_Conc ~ Day_numeric, data = Ecoli_select1_noNC)

# Generate 4 diagnostic plots (Residuals vs Fitted is first)
plot(Ecolimodel) 
BIC(Ecolimodel)
# To see only the Residuals vs Fitted plot:
plot(Ecolimodel, which = 1)
ggsave("/local1/workdir1/tw488/CEA_RData_RScript/Challenge_study/Check_linearity/Ecoliresiduals.tiff", units="in", width=20, height=15, dpi=500, compression = 'lzw')

# qqplotplot:
plot(Ecolimodel, which = 2)
ggsave("/local1/workdir1/tw488/CEA_RData_RScript/Challenge_study/Check_linearity/Ecoliqq.tiff", units="in", width=20, height=15, dpi=500, compression = 'lzw')

#divide Ecoli_noNC in to organic and conventional
Ecoli_select1_noNC_O<-Ecoli_select1_noNC%>%filter(Pond_Type=="Organic")
Ecoli_select1_noNC_C<-Ecoli_select1_noNC%>%filter(Pond_Type=="Conventional")

#linear model
Ecoli_select1_noNC_C_lm <- lm(Log_Average_Conc ~ Day_numeric, data = Ecoli_select1_noNC_C)
BIC(Ecoli_select1_noNC_C_lm)
Ecoli_select1_noNC_O_lm <- lm(Log_Average_Conc ~ Day_numeric, data = Ecoli_select1_noNC_O)
BIC(Ecoli_select1_noNC_O_lm)

#mafart conventional
mafart_weibull_Ecoli_C <- nls(Log_Average_Conc ~ log10_N0 - (Day_numeric/delta)^p, 
                   data =Ecoli_select1_noNC_C, 
                   start = list(log10_N0 = 4.612784, delta = 2, p = 1))


overview(mafart_weibull_Ecoli_C )
plotfit(mafart_weibull_Ecoli_C , smooth = TRUE)
BIC(mafart_weibull_Ecoli_C )
mafartparameters_Ecoli_C  <- summary(mafart_weibull_Ecoli_C )$coefficients
print(mafartparameters_Ecoli_C )
plot(fitted(mafart_weibull_Ecoli_C ), residuals(mafart_weibull_Ecoli_C ))
abline(h = 0, lty = 2, col = "red")

#mafart organic
mafart_weibull_Ecoli_O <- nls(Log_Average_Conc ~ log10_N0 - (Day_numeric/delta)^p, 
                              data =Ecoli_select1_noNC_O, 
                              start = list(log10_N0 = 4.612784, delta = 2, p = 1))


overview(mafart_weibull_Ecoli_O)
plotfit(mafart_weibull_Ecoli_O, smooth = TRUE)
BIC(mafart_weibull_Ecoli_O)
mafartparameters_Ecoli_O  <- summary(mafart_weibull_Ecoli_O)$coefficients
print(mafartparameters_Ecoli_O )
plot(fitted(mafart_weibull_Ecoli_O), residuals(mafart_weibull_Ecoli_O))
abline(h = 0, lty = 2, col = "red")

#get overall averages by organic/conventional status
Ecoli_bypond <- Ecoli_select1 %>% 
  filter(Aliquot!="NC")%>% 
  group_by(Day,Pond_Type) %>% 
  mutate(Average_Log_Conc_byPond=(mean(Log_Average_Conc)))%>%
  mutate(log_SD_Average_Conc_byPond=sd(Log_Average_Conc)) %>%
  ungroup

Ecoli_bypond <- Ecoli_bypond %>%distinct(Day,Pond_Type, .keep_all = TRUE)

Ecoli_bypond_fortable <-Ecoli_bypond%>%select(Sample_ID, Bacteria,Average_Log_Conc_byPond,log_SD_Average_Conc_byPond)%>%
  mutate(Average_Log_Conc_byPond = round(Average_Log_Conc_byPond, digits = 2))%>%
  mutate(log_SD_Average_Conc_byPond = round(log_SD_Average_Conc_byPond, digits = 2))


#find average of all EHEC NC samples
Ecoli_select1_NC_overall_avg<-Ecoli_select1_NC%>%mutate(NC_overall_avg=log10(mean(Average_Conc)))


#combine all averaged quantities by pond
Ecoli_bypond_fortable 
Salm_A_C_bypond_fortable 
Lmono_A_C_bypond_fortable 
stacked_df <- bind_rows(Lmono_A_C_bypond_fortable, Salm_A_C_bypond_fortable )

stacked_df <- bind_rows(stacked_df, Ecoli_bypond_fortable )
write.csv(stacked_df, "/local1/workdir1/tw488/CEA_RData_RScript/Challenge_study/PLate_Count_Values_for_Paper/All_Pathogen_Log_Conc_by_Day.csv", row.names = FALSE)

## ---------------------------
## Allequash Zooplankton Competitive Release/ Environmental Correlations
##      *Goal to seee if more variation is explained from environemnt or from Daph. pulicaria.
##
## Author: Kelly A. Loria
## Date Created: 2020-06-22
## Email: kelly.loria@colorado.edu
##
## ---------------------------
## Load packages:
library(ggplot2)
library(lubridate)
library(tidyverse)
library(vegan)
library(dplyr)
library(lme4)
library(lmerTest)
#library(plyr)

## ---------------------------
# File path setup:
if (dir.exists('/Users/kellyloria/Documents/Johnson\ Lab/Allequash/2019DaphniaProject')){
  inputDir<- '/Users/kellyloria/Documents/Johnson\ Lab/Allequash/2019DaphniaProject/Daph_Data'
  outputDir<- '/Users/kellyloria/Documents/Johnson\ Lab/Allequash/2019DaphniaProject/Daph_Data/Rproject_Output' 
}

## ---------------------------
# I. Read in data/fix timestamp
df_ZCR <- read.csv(paste0(inputDir, "/Daphnia_within_year_MASTER.csv"), header=T)
summary(df_ZCR)

# Convert df.env$date dateTime string to R date structure (date)                                
tmpDateFormat<-"%m/%d/%y"
tmp1break_date<-as.Date(df_ZCR$Date,format=tmpDateFormat)
# Keep the new dates only if they all converted correctly
if(length(tmp1break_date) == length(tmp1break_date[!is.na(tmp1break_date)])){df_ZCR$Date <- tmp1break_date } else {print("Date conversion failed for dt1$break_date. Please inspect the data and do the date conversion yourself.")}                                                                    
rm(tmpDateFormat,tmp1break_date) 

# Subset for columns of interest:
names(df_ZCR)
df_ZCR1 <- df_ZCR %>% 
  subset(  
    select=c(Year, Date, DayOfYear, Total_N, Total_Adult_N, Adult_Uninfected, Adult_Infected,
             Prevalence, Uninfected_Ungravid, Uninfected_Gravid, Uninfected_Male, Uninfected_Ephippia,
             Infected_Ungravid, Infected_Gravid, Infected_Male, Infected_Ephippia, Neonates, DAGA, 
             CERI, BOS, HOLO, DIA, Infected_DAGA, Infected_HOLO, GravidN, EggRatio, InfectedN, 
             EggsPerInfected, TotalEggEstimate, AverageEggs, MeasuredUnN, UnBodySize..mm., 
             MeasuredInfN, InBodySize..mm., Secchi, DailyAirTemp, DailyWindSpeed, ZoopRichness, 
             ZoopDensity, BosminaDens, CerioDens, DmendoateDens, DpulicariaDens, DiaphanosomaDens, 
             HolopediumDens, MeasurementDepth, SurfacePH, HypolimnionPH, SurfaceDOC, HypolimnionDOC,
             SurfaceTOC, HypolimnionTOC, SurfaceTDN, HypolimnionTDN, SurfaceTDP, HypolimnionTDP,
             ChlorophyllMax, SurfaceTemp, HypolimnionTemp
    )) 
summary(df_ZCR1) 


# Look for Chl-a correlations with all zoop groups:
pChla_B <- ggplot(df_ZCR1, aes(x=ChlorophyllMax, y=BosminaDens, alpha = 0.5)) +
  geom_point()+ theme_bw() #+ geom_smooth(method="lm", se=F)

lmChla_B <- glm(BosminaDens~ChlorophyllMax, data = df_ZCR1) 
summary(lmChla_B) # not sig 

pChla_CD <- ggplot(df_ZCR1, aes(x=ChlorophyllMax, y=CerioDens, alpha = 0.5)) +
  geom_point()+ theme_bw() + geom_smooth(method="lm", se=F)

lmChla_CD<- glm(CerioDens~ChlorophyllMax, data = df_ZCR1) 
summary(lmChla_CD) # marginal sig 

pChla_DM <- ggplot(df_ZCR1, aes(x=ChlorophyllMax, y=DmendoateDens, alpha = 0.5)) +
  geom_point()+ theme_bw() 

lmChla_DM <- glm(DmendoateDens ~ ChlorophyllMax, data = df_ZCR1) 
summary(lmChla_DM) # not sig 

pChla_DPHAN <- ggplot(df_ZCR1, aes(x=ChlorophyllMax, y=DiaphanosomaDens, alpha = 0.5)) +
  geom_point()+ theme_bw()

lmChla_DPHAN <- glm(DiaphanosomaDens ~ ChlorophyllMax, data = df_ZCR1) 
summary(lmChla_DPHAN) # not sig

pChla_HG <- ggplot(df_ZCR1, aes(x=ChlorophyllMax, y=HolopediumDens, alpha = 0.5)) +
  geom_point()+ theme_bw()

lmChla_HG <- glm(HolopediumDens ~ ChlorophyllMax, data = df_ZCR1) 
summary(lmChla_HG) #not sig

pChla_DP <- ggplot(df_ZCR1, aes(x=ChlorophyllMax, y=DpulicariaDens, alpha = 0.5)) +
  geom_point()+ theme_bw() + geom_smooth(method="lm", se=F)

lmChla_DP <- glm(DpulicariaDens ~ ChlorophyllMax, data = df_ZCR1) 
summary(lmChla_DP) # sig

# Notes:
# Based off just chl-a data D.pulicaria may be most effective grazers. 
# Also positive relationship of Cerio on Chla, but no relationship btwn Cerio and D.pulicaria 

# Look for Water Temp correlations with all zoop groups/ incase of winter dynamics:
pWT_B <- ggplot(df_ZCR1, aes(x=HypolimnionTemp, y=BosminaDens, alpha = 0.5)) +
  geom_point()+ theme_bw() 

lmWT_B <- glm(BosminaDens~HypolimnionTemp, data = df_ZCR1) 
summary(lmWT_B) # not sig 

pWT_CD <- ggplot(df_ZCR1, aes(x=HypolimnionTemp, y=CerioDens, alpha = 0.5)) +
  geom_point()+ theme_bw() + geom_smooth(method="lm", se=F)

lmWT_CD<- glm(CerioDens~HypolimnionTemp, data = df_ZCR1) 
summary(lmWT_CD) # Sig +

pWT_DM <- ggplot(df_ZCR1, aes(x=HypolimnionTemp, y=DmendoateDens, alpha = 0.5)) +
  geom_point()+ theme_bw() 

lmWT_DM <- glm(DmendoateDens ~ HypolimnionTemp, data = df_ZCR1) 
summary(lmWT_DM) # not sig 

pWT_DPHAN <- ggplot(df_ZCR1, aes(x=HypolimnionTemp, y=DiaphanosomaDens, alpha = 0.5)) +
  geom_point()+ theme_bw() + geom_smooth(method="lm", se=F)

lmWT_DPHAN <- glm(DiaphanosomaDens ~ HypolimnionTemp, data = df_ZCR1) 
summary(lmWT_DPHAN) # sig +

pWT_HG <- ggplot(df_ZCR1, aes(x=HypolimnionTemp, y=HolopediumDens, alpha = 0.5)) +
  geom_point()+ theme_bw()

lmWT_HG <- glm(HolopediumDens ~ HypolimnionTemp, data = df_ZCR1) 
summary(lmWT_HG) #not sig

pWT_DP <- ggplot(df_ZCR1, aes(x=HypolimnionTemp, y=DpulicariaDens, alpha = 0.5)) +
  geom_point()+ theme_bw() + geom_smooth(method="lm", se=F)

lmWT_DP <- glm(DpulicariaDens ~ HypolimnionTemp, data = df_ZCR1) 
summary(lmWT_DP) # sig -

# Notes:
# Based off hypo water temperature, cerio and Diaphanosoma are positively associated
#   while D.pulicaria is negatively associated

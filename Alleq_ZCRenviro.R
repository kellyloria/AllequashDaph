## ---------------------------
## Allequash Zooplankton Competitive Release/ Environmental Correlations
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
pDP_B <- ggplot(df_ZCR1, aes(x=DpulicariaDens, y=BosminaDens, alpha = 0.5)) +
  geom_point()+ theme_bw() + geom_smooth(method="lm", se=F)

lmDP_B <- glm(BosminaDens~DpulicariaDens, data = df_ZCR1) 
summary(lmDP_B) # marginal sig 

pDP_CD <- ggplot(df_ZCR1, aes(x=DpulicariaDens, y=CerioDens, alpha = 0.5)) +
  geom_point()+ theme_bw() 

lmDP_CD<- glm(CerioDens~DpulicariaDens, data = df_ZCR1) 
summary(lmDP_CD) # not sig 

pDP_DM <- ggplot(df_ZCR1, aes(x=DpulicariaDens, y=DmendoateDens, alpha = 0.5)) +
  geom_point()+ theme_bw() 

lmDP_DM <- glm(DmendoateDens ~ DpulicariaDens, data = df_ZCR1) 
summary(lmDP_DM) # not sig 

pDP_DPHAN <- ggplot(df_ZCR1, aes(x=DpulicariaDens, y=DiaphanosomaDens, alpha = 0.5)) +
  geom_point()+ theme_bw()

lmDP_DPHAN <- glm(DiaphanosomaDens ~ DpulicariaDens, data = df_ZCR1) 
summary(lmDP_DPHAN) #not sig

pDP_HG <- ggplot(df_ZCR1, aes(x=DpulicariaDens, y=HolopediumDens, alpha = 0.5)) +
  geom_point()+ theme_bw()

lmDP_HG <- glm(HolopediumDens ~ DpulicariaDens, data = df_ZCR1) 
summary(lmDP_HG) #not sig



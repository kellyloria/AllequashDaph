## ---------------------------
## Allequash Zooplankton Competitive Release
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

## ---------------------------
# II. Data exploration:
# Start with density overtime to see if there are any obvious visual tradeoffs. 

# Basic figure code/ convert from wide to long:
df_ZCRlong <- df_ZCR1 %>% 
  subset(Year>1999 & Year<2017, #removed 1999 and 2017 due to missing data 
    select=c(Year, Date, DayOfYear, BosminaDens, CerioDens, DmendoateDens, 
             DpulicariaDens, DiaphanosomaDens, HolopediumDens 
    )) %>%
  gather(taxon, density, BosminaDens:HolopediumDens, factor_key=TRUE)

# Plot of long format data
pZCR1_den <- ggplot(df_ZCRlong, aes(x=DayOfYear, y=log10(density+1), # log10 just to be able to see the data 
                                    shape = taxon, colour =taxon, alpha = 0.5)) +
  geom_point()+ theme_bw() + xlab("Days post ice-off") + facet_grid(Year~.)
#ggsave(paste0(outputDir,("/Alleq_pZCR1_den.pdf")), pZCR1_den, scale = 1.5, width = 10, height = 30, units = c("cm"), dpi = 500)
# looks like there might be some trade off-- but it could just be due to ice cover 

# Look at DpulicariaDens vs other groups 
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

# Issues: 
# 1. Might not be enough data positive obs of taxa to use ordination 

## ---------------------------
# II. Data exploration:
#  Modelling bimass of certain groups with envir




## ---------------------------
# EXTRA SCRAP CODE

# From DF data:
pZCR1_den2 <- ggplot(df_ZCR1) +
  geom_point(aes(x=DayOfYear, y=BosminaDens), shape =19, colour ="#7ba3b3", alpha = 0.5) +
  geom_point(aes(x=DayOfYear, y=CerioDens), shape =19, colour ="#7bb398", alpha = 0.5) +
  geom_point(aes(x=DayOfYear, y=DmendoateDens), shape =19, colour ="#95b37b", alpha = 0.5) +
  geom_point(aes(x=DayOfYear, y=DpulicariaDens), shape =17, colour ="#b81500", alpha = 0.5) +
  geom_point(aes(x=DayOfYear, y=DiaphanosomaDens), shape =19, colour ="#7b8ab3", alpha = 0.5) +
  geom_point(aes(x=DayOfYear, y=HolopediumDens), shape =19, colour ="#7f7bb3", alpha = 0.5) +
  theme_classic() + xlab("Days post ice-off") #+ facet_grid(Year~.)




## ---------------------------
## Allequash Zooplankton Environmental Effects
##      *Time series patterns for water quality 
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


# Basic figure code/ convert from wide to long:
df_Physlong <- df_ZCR1 %>% 
  subset(#Year>1999 & Year<2017, #removed 1999 and 2017 due to missing data 
         select=c(Year, Date, DayOfYear, Secchi, DailyAirTemp, DailyWindSpeed, 
                  SurfacePH, HypolimnionPH, SurfaceDOC, HypolimnionDOC, 
                  SurfaceTDN, HypolimnionTDN, SurfaceTDP, HypolimnionTDP,
                  ChlorophyllMax, SurfaceTemp, HypolimnionTemp
         )) %>%
  gather(Envir.sp, value, Secchi:HypolimnionTemp, factor_key=TRUE)


# Plot of long format data
pPhys_TS <- ggplot(df_Physlong, aes(x=DayOfYear, y=log10(value+1), # log10 just to be able to see the data 
                                    shape = Envir.sp, colour =Envir.sp, alpha = 0.5)) +
  geom_point()+ theme_bw() + xlab("Days post ice-off") + facet_grid(Year~.) ## NOT very helpful



## ---------------------------
# PCA to figure out what level of variable explains most variation

df_Phys <- df_ZCR1 %>% 
  subset(#Year>1999 & Year<2017, #removed 1999 and 2017 due to missing data 
    select=c(Year, Date, DayOfYear, Secchi, DailyAirTemp, DailyWindSpeed, 
             SurfacePH, HypolimnionPH, SurfaceDOC, HypolimnionDOC, 
             SurfaceTDN, HypolimnionTDN, SurfaceTDP, HypolimnionTDP,
             ChlorophyllMax, SurfaceTemp, HypolimnionTemp
    ))
# I feel like there might be some missing data 
###### %$#^%%$%
# Stoping here to see if some of Agg code has better resolution maybe more overlap?


all.dat.maxns1 %>%
  metaMDS(trace = T) %>%
  ordiplot(type = "none") %>%
  text("sites")

PCA <- rda(all.dat.maxns1, scale = FALSE)
# Use scale = TRUE if your variables are on different scales (e.g. for abiotic variables).
# Here, all species are measured on the same scale 
# So use scale = FALSE

# Now plot a bar plot of relative eigenvalues. This is the percentage variance explained by each axis
barplot(as.vector(PCA$CA$eig)/sum(PCA$CA$eig)) 
# How much of the variance in our dataset is explained by the first principal component?

# Calculate the percent of variance explained by first two axes
sum((as.vector(PCA$CA$eig)/sum(PCA$CA$eig))[1:2]) # 0.977%

# Now, we`ll plot our results with the plot function
plot(PCA)
plot(PCA, display = "sites", type = "points")
plot(PCA, display = "species", type = "text")


#   2. Start with average:
all.dat.av.ny2 %>%
  metaMDS(trace = T) %>%
  ordiplot(type = "none") %>%
  text("sites")

PCA <- rda(all.dat.av.ny2, scale = FALSE)
# Use scale = TRUE if your variables are on different scales (e.g. for abiotic variables).
# Here, all species are measured on the same scale 
# So use scale = FALSE

# Now plot a bar plot of relative eigenvalues. This is the percentage variance explained by each axis
barplot(as.vector(PCA$CA$eig)/sum(PCA$CA$eig)) 
# How much of the variance in our dataset is explained by the first principal component?

# Calculate the percent of variance explained by first two axes
sum((as.vector(PCA$CA$eig)/sum(PCA$CA$eig))[1:2]) # 0.967982%

# Now, we`ll plot our results with the plot function
plot(PCA)
plot(PCA, display = "sites", type = "points")
plot(PCA, display = "species", type = "text")

# Recap TDN big driver regardless of how data is spliced. Without TDN and TDP, O2sat and chlora standout


## ---------------------------
## Aggregation of NTL Allequash data
##
## Author: Kelly A. Loria
## Date Created: 2019-03-02 & updated 2020-03-19
## Email: kelly.loria@colorado.edu
##
## ---------------------------
## Load packages:
library(ggplot2)
library(dplyr)
library(tidyverse)

setwd("~/Documents/Johnson Lab/Allequash")

## ---------------------------
# I. Grab and clean data:
#   1. Start with Secchi Depth: 
AL_sechi <- read.csv("ntl31_v6_secchi.csv")
summary(AL_sechi)

#	secview	= secchi depth with viewer
# secnview = secchi depth without viewer

AL_sechi$date <- as.Date.character(AL_sechi$sampledate, format="%m/%d/%y")
range(AL_sechi$date)
summary(AL_secchi.s)

AL_secchi.s <- subset(AL_sechi, select=c(year4, date, daynum, secview,
                                         airtemp, windspd))
colnames(AL_secchi.s)[1] = "year"
colnames(AL_secchi.s)[2] = "date"
colnames(AL_secchi.s)[3] = "day.of.year"
colnames(AL_secchi.s)[4] = "secchi"
colnames(AL_secchi.s)[5] = "daily.airtemp"
colnames(AL_secchi.s)[6] = "daily.windsp"

# New .csv with relevant column headers: 
# write.csv(AL_secchi.s, "NTL.Alleq_secchi.csv")

## ---------------------------
#   2. Pull in water quality parameters
#     * Data source: https://portal.edirepository.org/nis/metadataviewer?packageid=knb-lter-ntl.1.49
#     * There are replicates of data per dates
AL_WQ <- read.csv("chemphys.csv")
summary(AL_WQ)

AL_WQ$date <- as.Date.character(AL_WQ$sampledate, format="%m/%d/%y")
range(AL_WQ$date)
summary(AL_WQ.s)

AL_WQ.s <- subset(AL_WQ, select=c(year4, date, daynum, depth,
                                  ph, doc, toc, totnf, totpf))
colnames(AL_WQ.s)[1] = "year"
colnames(AL_WQ.s)[2] = "date"
colnames(AL_WQ.s)[3] = "day.of.year"
colnames(AL_WQ.s)[8] = "TDN"
colnames(AL_WQ.s)[9] = "TDP"

# some values likely not real and left over from flagging like -99.99
AL_WQ.s1 <- AL_WQ.s %>% 
  subset(doc >= 0 & toc >= 0 & toc <= 99)
summary(AL_WQ.s1)


AL_WQ.s2 <- AL_WQ.s1 %>% 
  subset(TDN >= 0 & TDP >= 0 )
summary(AL_WQ.s2)

# New .csv with relevant column headers: 
# write.csv(AL_WQ.s, "NTL.Alleq_OrgNut.csv")

## add OrgNuts to secchi ##
AL_WQ1 <- left_join(AL_WQ.s2, AL_secchi.s[c("date", "secchi", "daily.airtemp", "daily.windsp")],
                    by = c("date" = "date"))
summary(AL_WQ1)

## ---------------------------
#   3. Pull in chl-a 
AL_chla <- read.csv("chlorophyll_nl_q1.csv")
summary(AL_chla)

AL_chla$date <- as.Date.character(AL_chla$sampledate, format="%m/%d/%y")
range(AL_chla$date)

AL_chla$Julian <- as.POSIXlt( as.Date(AL_chla$date,format= "%m/%d/%y"), 
                              format= "%Y-%m-%d")$yday

# unclear values high values
# https://www.rmbel.info/primer/chlorophyll-a/
# 150 is hyper eutrophic so restricting for values greater than 200
AL_chla.Q <- AL_chla %>% 
  subset(chlor >= 0 & chlor <= 190 & phaeo >=0 & phaeo <=190)
summary(AL_chla.Q)

AL_WQ2 <- left_join(AL_WQ1, AL_chla.Q[c("date", "depth", "chlor", "phaeo")],
                    by = c("date" = "date", "depth" = "depth"))
summary(AL_WQ2)

AL.secchi.p <- ggplot(AL_WQ2) +
  geom_point(aes(x=day.of.year, y=((secchi))), shape =19, colour ="#ad8800", alpha = 0.5) +
  geom_line(aes(x=day.of.year, y=((secchi))), colour ="#ad8800", alpha = 0.5) +
  theme_classic() + ylab("Secchi Depth (meters)") +xlab("Date") +
  facet_wrap(~year)

AL.chla.p <- ggplot(AL_WQ2) +
  geom_point(aes(x=day.of.year, y=((chlor))), shape =19, colour ="#ad8800", alpha = 0.5) +
  geom_line(aes(x=day.of.year, y=((chlor))), colour ="#ad8800", alpha = 0.5) +
  theme_classic() + ylab("Chla (ugL)") +xlab("Date") +
  facet_wrap(~year) # values over 200 might not be real

# write.csv(AL_WQ2, "NTL.Alleq_OrgNutChla.csv")

## ---------------------------
#   4. Also add in DO + water temp
AL_WT <- read.csv("chemphys.csv")
summary(AL_WT)

AL_WT$date <- as.Date.character(AL_WT$sampledate, format="%Y-%m-%d")
range(AL_WT$date)
summary(AL_WT)

AL_WT.s <- subset(AL_WT, select=c(year4, date, daynum, depth,
                                  wtemp, o2, o2sat, light))

# some values likely not real and left over from flagging like -99.99
AL_WQ.s1 <- AL_WQ.s %>% 
  subset(doc >= 0 & toc >= 0 & toc <= 99)
summary(AL_WQ.s1)

AL.chla.p <- ggplot(AL_WT) +
  geom_point(aes(x=daynum, y=((wtemp))), shape =19, colour ="#ad8800", alpha = 0.5) +
  geom_line(aes(x=daynum, y=((wtemp))), colour ="#ad8800", alpha = 0.5) +
  theme_classic() + ylab("WaterT (C)") +xlab("Date") +
  facet_wrap(~year4) # values over 200 might not be real

AL.chla.p <- ggplot(AL_WT) +
  geom_point(aes(x=daynum, y=((o2))), shape =19, colour ="#ad8800", alpha = 0.5) +
  geom_line(aes(x=daynum, y=((o2))), colour ="#ad8800", alpha = 0.5) +
  theme_classic() + ylab("WaterT (C)") +xlab("Date") +
  facet_wrap(~year4) # values over 200 might not be real

AL_WQ3 <- left_join(AL_WQ2, AL_WT.s[c("date", "depth", "wtemp", "o2", "o2sat", "light")],
                    by = c("date" = "date", "depth" = "depth"))
summary(AL_WQ3)

# write.csv(AL_WQ3, "NTL.Alleq_OrgNutChlaTemp.csv")

## ---------------------------
# II. 2020-03-03 timeseries model 
#     * Double check the lake physical models
#     * Run auto regressive models on data set - most auto regressive models can't handle 
library(glmmTMB)
library(dplyr) # need to group by year so we only implement by year

#data <- read.csv(file.choose())
summary(data)

#phy_data <- read.csv(file.choose())
summary(phy_data)

data$date1 <- as.Date(data$Date, format="%m/%d/%y")
range(data$date1)

#add in raw phys
data2 <- left_join(data, AL_secchi.s[c("date","secview",
                                       "airtemp", "windspd")],
                   by = c("date1" = "date"))

write.csv(data2, "zoop2000prac_secchi.csv")

mdl <- glmmTMB(response ~ scale(predictor) + ar1(time + 0 | grpvar), family = binomial, ziformula=~0, data = data)

lag1 <- data %>%
  group_by(Year) %>%
  mutate(lagN1 = lag(infected, 1, order_by=julianall))

## ---------------------------
## III. Calculate average and max per depth/day for 
##      PCA to figure out what level of variable explains most variation
library(naniar)
library(vegan)
#   1. Read in past year's data - here 2018 summer
all.dat <- read.csv("NTL.Alleq_OrgNutChlaTemp.csv", header=T)

#   2. Fix timestamp - so it is no longer a character:
all.dat$date1 <- as.Date(all.dat$date, format="%m/%d/%y ")
range(all.dat$date1)
names(all.dat)

#   3. Calculate max per depth/day:
all.dat.maxny <- all.dat %>% 
  group_by(date1) %>% 
  summarise("ph.m" = max(ph, na.rm = T),
            "doc.m"= max(doc, na.rm = T),
            "toc.m"= max(toc, na.rm = T),
            #"TDN.m"= max(TDN, na.rm = T),
            #"TDP.m"= max(TDP, na.rm = T),
            "secchi.m"= max(secchi, na.rm = T),
            "daily.airtemp.m"= max(daily.airtemp, na.rm = T),
            "daily.windsp.m"= max(daily.windsp, na.rm = T),
            "chlor.m"= max(chlor, na.rm = T),
            "wtemp.m"= max(wtemp, na.rm = T),
            "o2.m"= max(o2, na.rm = T),
            "o2sat.m"= max(o2sat, na.rm = T),
            "light.m"= max(light, na.rm = T))
summary(all.dat.maxny)

# eliminate all -inf
all.dat.maxny1 <- all.dat.maxns[!is.na(x) & !is.infinite(x)]

#   4. Select only numeric values:
all.dat.maxns <- select(all.dat.maxny, 2:11)
all.dat.maxns1 <- all.dat.maxns %>% 
  filter_all(all_vars(is.finite(.)))

#   5. Calculate average per depth/day
all.dat.av.ny <- all.dat %>% 
  group_by(date1) %>% 
  summarise("ph.a" = mean(ph, na.rm = T),
            "doc.a"= mean(doc, na.rm = T),
            "toc.a"= mean(toc, na.rm = T),
            "TDN.a"= mean(TDN, na.rm = T),
            "TDP.a"= mean(TDP, na.rm = T),
            "secchi.a"= mean(secchi, na.rm = T),
            "daily.airtemp.a"= mean(daily.airtemp, na.rm = T),
            "daily.windsp.a"= mean(daily.windsp, na.rm = T),
            "chlor.a"= mean(chlor, na.rm = T),
            "wtemp.a"= mean(wtemp, na.rm = T),
            "o2.a"= mean(o2, na.rm = T),
            "o2sat.a"= mean(o2sat, na.rm = T),
            "light.a"= mean(light, na.rm = T))
summary(all.dat.av.ny)

all.dat.av.ny1 <- na.omit(all.dat.av.ny)

# eliminate all -inf
all.dat.maxny1 <- all.dat.maxns[!is.na(x) & !is.infinite(x)]

#   4. Select only numeric values:
all.dat.av.ny2 <- select(all.dat.av.ny1, 2:13)


## ---------------------------
## IV. PCA to figure out what level of variable explains most variation

#   1. Start with max:
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







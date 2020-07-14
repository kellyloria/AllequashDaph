## ---------------------------
## Allequash Ice thickness cleaning script
## ---------------------------
## Load packages:
library(ggplot2)
library(lubridate)
library(tidyverse)
library(vegan)
library(dplyr)
library(lme4)
library(lmerTest)

## ---------------------------
# File path setup:
if (dir.exists('/Users/kellyloria/Documents/Johnson\ Lab/Allequash/2019DaphniaProject')){
  inputDir<- '/Users/kellyloria/Documents/Johnson\ Lab/Allequash/2019DaphniaProject/Daph_Data'
  outputDir<- '/Users/kellyloria/Documents/Johnson\ Lab/Allequash/2019DaphniaProject/Daph_Data/Rproject_Output' 
}

## ---------------------------
# Ice thickness data: https://portal.edirepository.org/nis/codeGeneration?packageId=knb-lter-ntl.34.31&statisticalFileType=r

# Package ID: knb-lter-ntl.34.31 Cataloging System:https://pasta.edirepository.org.
# Data set title: North Temperate Lakes LTER: Snow and Ice Depth 1982 - current.
# Data set creator:  NTL Lead PI - University of Wisconsin 
# Data set creator:  John Magnuson - University of Wisconsin 
# Data set creator:  Stephen Carpenter - University of Wisconsin 
# Data set creator:  Emily Stanley - University of Wisconsin 
# Metadata Provider:  NTL Information Manager - University of Wisconsin 
# Contact:  NTL Information Manager -  University of Wisconsin  - ntl.infomgr@gmail.com
# Contact:  NTL Lead PI -  University of Wisconsin  - ntl.leadpi@gmail.com
# Stylesheet v2.11 for metadata conversion into program: John H. Porter, Univ. Virginia, jporter@virginia.edu 

inUrl1  <- "https://pasta.lternet.edu/package/data/eml/knb-lter-ntl/34/31/9af7f7d823fd8be3e4e31ccd7d4bb003" 
infile1 <- tempfile()
try(download.file(inUrl1,infile1,method="curl"))
if (is.na(file.size(infile1))) download.file(inUrl1,infile1,method="auto")


dt1 <-read.csv(infile1,header=F 
               ,skip=1
               ,sep=","  
               ,quot='"' 
               , col.names=c(
                 "lakeid",     
                 "year4",     
                 "daynum",     
                 "sampledate",     
                 "sta",     
                 "nsnow",     
                 "avsnow",     
                 "sdsnow",     
                 "wlevel",     
                 "totice",     
                 "nice",     
                 "whiteice",     
                 "blueice"    ), check.names=TRUE)

unlink(infile1)

# Fix any interval or ratio columns mistakenly read in as nominal and nominal columns read as numeric or dates read as strings

if (class(dt1$lakeid)!="factor") dt1$lakeid<- as.factor(dt1$lakeid)
if (class(dt1$year4)=="factor") dt1$year4 <-as.numeric(levels(dt1$year4))[as.integer(dt1$year4) ]               
if (class(dt1$year4)=="character") dt1$year4 <-as.numeric(dt1$year4)
if (class(dt1$daynum)=="factor") dt1$daynum <-as.numeric(levels(dt1$daynum))[as.integer(dt1$daynum) ]               
if (class(dt1$daynum)=="character") dt1$daynum <-as.numeric(dt1$daynum)                                   
# attempting to convert dt1$sampledate dateTime string to R date structure (date or POSIXct)                                
tmpDateFormat<-"%Y-%m-%d"
tmp1sampledate<-as.Date(dt1$sampledate,format=tmpDateFormat)
# Keep the new dates only if they all converted correctly
if(length(tmp1sampledate) == length(tmp1sampledate[!is.na(tmp1sampledate)])){dt1$sampledate <- tmp1sampledate } else {print("Date conversion failed for dt1$sampledate. Please inspect the data and do the date conversion yourself.")}                                                                    
rm(tmpDateFormat,tmp1sampledate) 
if (class(dt1$sta)!="factor") dt1$sta<- as.factor(dt1$sta)
if (class(dt1$nsnow)=="factor") dt1$nsnow <-as.numeric(levels(dt1$nsnow))[as.integer(dt1$nsnow) ]               
if (class(dt1$nsnow)=="character") dt1$nsnow <-as.numeric(dt1$nsnow)
if (class(dt1$avsnow)=="factor") dt1$avsnow <-as.numeric(levels(dt1$avsnow))[as.integer(dt1$avsnow) ]               
if (class(dt1$avsnow)=="character") dt1$avsnow <-as.numeric(dt1$avsnow)
if (class(dt1$sdsnow)=="factor") dt1$sdsnow <-as.numeric(levels(dt1$sdsnow))[as.integer(dt1$sdsnow) ]               
if (class(dt1$sdsnow)=="character") dt1$sdsnow <-as.numeric(dt1$sdsnow)
if (class(dt1$wlevel)=="factor") dt1$wlevel <-as.numeric(levels(dt1$wlevel))[as.integer(dt1$wlevel) ]               
if (class(dt1$wlevel)=="character") dt1$wlevel <-as.numeric(dt1$wlevel)
if (class(dt1$totice)=="factor") dt1$totice <-as.numeric(levels(dt1$totice))[as.integer(dt1$totice) ]               
if (class(dt1$totice)=="character") dt1$totice <-as.numeric(dt1$totice)
if (class(dt1$nice)=="factor") dt1$nice <-as.numeric(levels(dt1$nice))[as.integer(dt1$nice) ]               
if (class(dt1$nice)=="character") dt1$nice <-as.numeric(dt1$nice)
if (class(dt1$whiteice)=="factor") dt1$whiteice <-as.numeric(levels(dt1$whiteice))[as.integer(dt1$whiteice) ]               
if (class(dt1$whiteice)=="character") dt1$whiteice <-as.numeric(dt1$whiteice)
if (class(dt1$blueice)=="factor") dt1$blueice <-as.numeric(levels(dt1$blueice))[as.integer(dt1$blueice) ]               
if (class(dt1$blueice)=="character") dt1$blueice <-as.numeric(dt1$blueice)

# Here is the structure of the input data frame:
str(dt1)                            
attach(dt1)                            
# The analyses below are basic descriptions of the variables. After testing, they should be replaced.                 
summary(as.factor(dt1$lakeid)) 

dt3 <- dt1 %>% 
  subset(lakeid == "AL")
summary(dt3)

# Calculate mean ice-thickness
dt3_summary = dt3 %>% 
  group_by(year4)%>% #this will get the nearest 15, but could be fewer if some are missing OR >35C, I think (?) the 35 are bogus so that is ok but you could
  mutate(mn_iceT=mean(totice),           # also filter out the NAs and >35s if you wanted to always have 15 values in your rolling window after removing bad values
         sd_iceT=sd(totice)) %>%
  mutate(mn_iceW=mean(whiteice),           # also filter out the NAs and >35s if you wanted to always have 15 values in your rolling window after removing bad values
         sd_iceW=sd(whiteice))%>%
  mutate(mn_iceB=mean(blueice),           # also filter out the NAs and >35s if you wanted to always have 15 values in your rolling window after removing bad values
         sd_iceB=sd(blueice))

## ---------------------------
###
# Package ID: knb-lter-ntl.32.27 Cataloging System:https://pasta.edirepository.org.
# Data set title: North Temperate Lakes LTER: Ice Duration - Trout Lake Area 1981 - current.
# Data set creator:  NTL Lead PI - University of Wisconsin 
# Data set creator:  North Temperate Lakes LTER -  
# Data set creator:  John Magnuson - University of Wisconsin 
# Data set creator:  Stephen Carpenter - University of Wisconsin 
# Data set creator:  Emily Stanley - University of Wisconsin 
# Metadata Provider:  NTL Information Manager - University of Wisconsin 
# Contact:  NTL Information Manager -  University of Wisconsin  - ntl.infomgr@gmail.com
# Contact:  NTL Lead PI -  University of Wisconsin  - ntl.leadpi@gmail.com
# Stylesheet v2.11 for metadata conversion into program: John H. Porter, Univ. Virginia, jporter@virginia.edu 

inUrl1  <- "https://pasta.lternet.edu/package/data/eml/knb-lter-ntl/32/27/e57a6b46a237355214844e2c76fa8aa5" 
infile1 <- tempfile()
try(download.file(inUrl1,infile1,method="curl"))
if (is.na(file.size(infile1))) download.file(inUrl1,infile1,method="auto")


dt1 <-read.csv(infile1,header=F 
               ,skip=1
               ,sep=","  
               ,quot='"' 
               , col.names=c(
                 "lakeid",     
                 "sta",     
                 "year",     
                 "lastice",     
                 "datelastice",     
                 "firstopen",     
                 "datefirstopen",     
                 "lastopen",     
                 "datelastopen",     
                 "firstice",     
                 "datefirstice"    ), check.names=TRUE)

unlink(infile1)

# Fix any interval or ratio columns mistakenly read in as nominal and nominal columns read as numeric or dates read as strings

if (class(dt1$lakeid)!="factor") dt1$lakeid<- as.factor(dt1$lakeid)
if (class(dt1$sta)!="factor") dt1$sta<- as.factor(dt1$sta)
if (class(dt1$year)=="factor") dt1$year <-as.numeric(levels(dt1$year))[as.integer(dt1$year) ]               
if (class(dt1$year)=="character") dt1$year <-as.numeric(dt1$year)
if (class(dt1$lastice)!="factor") dt1$lastice<- as.factor(dt1$lastice)                                   
# attempting to convert dt1$datelastice dateTime string to R date structure (date or POSIXct)                                
tmpDateFormat<-"%Y-%m-%d"
tmp1datelastice<-as.Date(dt1$datelastice,format=tmpDateFormat)
# Keep the new dates only if they all converted correctly
if(length(tmp1datelastice) == length(tmp1datelastice[!is.na(tmp1datelastice)])){dt1$datelastice <- tmp1datelastice } else {print("Date conversion failed for dt1$datelastice. Please inspect the data and do the date conversion yourself.")}                                                                    
rm(tmpDateFormat,tmp1datelastice) 
if (class(dt1$firstopen)!="factor") dt1$firstopen<- as.factor(dt1$firstopen)                                   
# attempting to convert dt1$datefirstopen dateTime string to R date structure (date or POSIXct)                                
tmpDateFormat<-"%Y-%m-%d"
tmp1datefirstopen<-as.Date(dt1$datefirstopen,format=tmpDateFormat)
# Keep the new dates only if they all converted correctly
if(length(tmp1datefirstopen) == length(tmp1datefirstopen[!is.na(tmp1datefirstopen)])){dt1$datefirstopen <- tmp1datefirstopen } else {print("Date conversion failed for dt1$datefirstopen. Please inspect the data and do the date conversion yourself.")}                                                                    
rm(tmpDateFormat,tmp1datefirstopen) 
if (class(dt1$lastopen)!="factor") dt1$lastopen<- as.factor(dt1$lastopen)                                   
# attempting to convert dt1$datelastopen dateTime string to R date structure (date or POSIXct)                                
tmpDateFormat<-"%Y-%m-%d"
tmp1datelastopen<-as.Date(dt1$datelastopen,format=tmpDateFormat)
# Keep the new dates only if they all converted correctly
if(length(tmp1datelastopen) == length(tmp1datelastopen[!is.na(tmp1datelastopen)])){dt1$datelastopen <- tmp1datelastopen } else {print("Date conversion failed for dt1$datelastopen. Please inspect the data and do the date conversion yourself.")}                                                                    
rm(tmpDateFormat,tmp1datelastopen) 
if (class(dt1$firstice)!="factor") dt1$firstice<- as.factor(dt1$firstice)                                   
# attempting to convert dt1$datefirstice dateTime string to R date structure (date or POSIXct)                                
tmpDateFormat<-"%Y-%m-%d"
tmp1datefirstice<-as.Date(dt1$datefirstice,format=tmpDateFormat)
# Keep the new dates only if they all converted correctly
if(length(tmp1datefirstice) == length(tmp1datefirstice[!is.na(tmp1datefirstice)])){dt1$datefirstice <- tmp1datefirstice } else {print("Date conversion failed for dt1$datefirstice. Please inspect the data and do the date conversion yourself.")}                                                                    
rm(tmpDateFormat,tmp1datefirstice) 

dt2 <- dt1 %>% 
  subset(lakeid == "AL")
summary(dt2)    


summary(dt3_summary)

dt4 <- left_join(dt2, dt3_summary[c("year4", "mn_iceT", "sd_iceT")],
                          by = c("year" = "year4"))

#   5. Double chec for duplicated values:
dt4%>%select(year)%>%duplicated()%>%sum() # 2 dups

View(dt4%>%
       inner_join(
         dt4 %>%
           group_by(year) %>%
           summarize(ct=dplyr::n())%>% filter(ct>1)))

# Remove values:
dt4 = dt4 %>%
  distinct(year, .keep_all = TRUE)

dt4%>%select(year)%>%duplicated()%>%sum() 


# write.csv(dt4, paste0(outputDir, "Alleq_Ice.csv")) 


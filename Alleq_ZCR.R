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
if (dir.exists('/Users/kellyloria/Documents/Johnson\ Lab/Allequash/2019\ Daphnia\ project/Daph_Data/Daphnia_within_year_MASTER.cs')){
  inputDir<- '/Users/kellyloria/Documents/Publications/GL4\ Overtime/'
  outputDir<- '/Users/kellyloria/Desktop/' 
}

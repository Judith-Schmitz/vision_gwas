
library(foreign) # for read.spss()
library(tidyverse) # for read.spss()
library(tidylog)
library(xlsx) # export dataframes to xlsx files
library(summarytools)
library(BBmisc) # for z standardisation
library(RNOmni) # for rank normalization
library(ggpubr)
library(nortest)
library(Hmisc)
library(corrplot)
library(RColorBrewer)
library(psych)

read.alspac <- function(mypath = 'data/Paracchini_03Apr19.sav') {
  
  # Purpose:  Import an Alspac dataset from SPSS and rename the columns to a clearer syntax 
  #           
  # Inputs:   path: .sav location in the file system 
  #
  # Output:   one dataframe with the Alspac data
  #
  
  default.opt <- options()
  options(warn = -1) # Avoid warnings "record type 7, subtype 22" of read.spss() function
  
  input.alspac <- read.spss(file = mypath,
                            use.value.labels = FALSE, # there are levels with undefined labels in the SPSS file
                            to.data.frame = TRUE,
                            use.missings = TRUE)
  
  # Remove labels
  for (var in colnames(input.alspac)) {
    attr(input.alspac[,deparse(as.name(var))], "value.labels") <- NULL
  }
  
  # Choose variables
  alspac <- input.alspac[,c('cid_417b', # unique pregnancy identifier (cid_417b / cidB823)
                            'qlet', # Identify children from multiple births
                            'kz021', # sex; 1 = male, 2 = female, -1 = not known
                            
                            # Age at 7 years
                            'f7003a', # age in days
                            'f7003b', # age in weeks
                            'f7003c', # age in months
                            
                            # Age at 11 years
                            'fe003a', # age in days
                            'fe003b', # age in weeks
                            'fe003c', # age in months
                            
                            # hearing
                            'f7hs017', # air conduction right
                            'f7hs027', # air conduction left
                            'f7hs035', # hearing impairment
                            'f7hs036', # hearing loss
                            'f7hs037', # high freq hearing loss
                            'f7hs071', # hearing letter given
                            
                            # vision
                            'fevs187', # best visual acuity
                            'fevs184', # difference
                            
                            # exclusion criteria controls
                            'ku991b', # age in weeks at 9 years
                            'kr803a', # ADHD diagnosis
                            'ku360', # ASD
                            'ku503b', # CCC 1 - autism
                            'ku504b', # CCC 2
                            'ku505b', # CCC 3
                            'ku506b', # CCC 4
                            'ku507b', # CCC 5
                            'ku508b', # CCC 6
                            'ku509b', # CCC 7
                            'f8sl105', # nonword rep - sli
                            'f8sl040', # WOLD comprehension - sli
                            'kr610', # language therapy - sli
                            
                            'sa036a', # Child has ever had sensory impairment (visual) (2=no)
                            
                            'c645a', # Mother highest educational qualification at pregnancy
                            'f8003b', # age in weeks F8
                            'f8ws112', # WISC total IQ (age corrected) - age 8
                            'f8ws110', # WISC verbal IQ (age corrected) - age 8
                            'f8ws111', # WISC performance IQ (age corrected) - age 8
                            'f8ws052', # WISC arithmetic (age corrected) - age 8
                            'f7ws076' # reading score at age 7
                            
  )]
  
  # Change names
  names(alspac) <- c("id", #cidB823
                     "id.twins", #qlet
                     "sex", #kz021
                     
                     "age.days.7", #f7003a  at 7
                     "age.weeks.7", #f7003b at 7
                     "age.months.7", #f7003c at 7 

                     "age.days.11",
                     "age.weeks.11",
                     "age.months.11", 
                     
                     "ear.right", 
                     "ear.left",
                     "hear.imp",
                     "hear.loss",
                     "hf.hear.loss",
                     "hear.letter",
                     
                     "best.eye",
                     "eye.diff",
                     
                     "age.weeks.9", #ku991b
                     "adhd", #kr803a
                     "asd", #ku360
                     "CCC1", #ku503b
                     "CCC2", #ku504b
                     "CCC3", #ku505b
                     "CCC4", #ku506b
                     "CCC5", #ku507b
                     "CCC6", #ku508b
                     "CCC7", #ku509b
                     "nonword.rep", #f8sl105
                     "wold.comp", #f8sl040
                     "speech.therapy", #kr610
                     
                     "sensory.impairment",
                     
                     "mother.education", #c645a
                     "age.weeks.8", #f8003b at 8
                     "wisc.total.iq", #f8ws112 at 8
                     "wisc.verbal.iq", #f8ws110 at 8
                     "wisc.performance.iq", #f8ws111 at 8
                     "wisc.arithmetic", #f8ws052 at 8
                     "read.7"
  )
  
  alspac <- alspac %>% mutate(
    ID_1 = as.character(paste0(id, id.twins)),
    sum.ccc = CCC1 + CCC2 + CCC3 + CCC4 + CCC5 + CCC6 + CCC7)

  options(default.opt) #Set back the default values
  
  return(alspac)
}


read.alspac.newfile <- function(mypath = 'data/Paracchini_12Feb20.sav') {
  
  # Purpose:  Import an Alspac dataset from SPSS and rename the columns to a clearer syntax 
  #           
  # Inputs:   path: .sav location in the file system 
  #
  # Output:   one dataframe with the Alspac data
  #
  
  default.opt <- options()
  options(warn = -1) # Avoid warnings "record type 7, subtype 22" of read.spss() function
  
  input.alspac <- read.spss(file = mypath,
                            use.value.labels = FALSE, # there are levels with undefined labels in the SPSS file
                            to.data.frame = TRUE,
                            use.missings = TRUE)
  
  # Remove labels
  for (var in colnames(input.alspac)) {
    attr(input.alspac[,deparse(as.name(var))], "value.labels") <- NULL
  }
  
  # Choose variables
  alspac <- input.alspac[,c('cid_417b', # unique pregnancy identifier (cid_417b / cidB823)
                            'qlet', # Identify children from multiple births
                            'ks4_ptscnewe' # Capped GCSE
                            
  )]
  
  # Change names
  names(alspac) <- c("id", #cidB823
                     "id.twins", #qlet
                     "GCSE" #ks4_ptscnewe
  )
  
  alspac <- alspac %>% mutate(
    ID_1 = as.character(paste0(id, id.twins))) %>%
    dplyr::select(
      ID_1,
      GCSE)
  
  alspac$GCSE[alspac$GCSE == "-10"] <- NA
  alspac$GCSE[alspac$GCSE == "-1"] <- NA
  
  options(default.opt) #Set back the default values
  
  return(alspac)
}


clean.alspac <- function(mydata) {
  
  # Purpose:  Filter subjects with sensory impairments and missing phenotype data
  #
  # Inputs:   mydata: an Alspac dataframe
  #
  # Output:   A modified Alspac dataframe
  #
  
  working.data <- mydata %>% 
    filter(!is.na(eye.diff) & !is.na(best.eye)) %>%
    filter(is.na(sensory.impairment) | sensory.impairment != 1 ) %>%
    filter(best.eye <= 0.3) %>%
    filter(eye.diff <= 0.2) %>%
    filter(!is.na(best.eye)) %>%
    # Mutate categorical variables as factors
    mutate(sex = as.factor(sex),
           id.twins = as.character(id.twins)) %>% 
    mutate_if(.predicate = is.factor, .funs = funs('levels<-'(., trimws(levels(.))))) #Trim factor-levels' string
  
  working.data[working.data == -9999] <- NA
  working.data[working.data == -69993] <- NA
  
  return(working.data)
  
}


normalize.alspac <- function(dataset) {

# Purpose:  Format Alspac measures
#
# Inputs:   dataset: cleaned Alspac dataframe 
#
# Output:   A normalized Alspac dataframe
#

  norm.data <- dataset %>% mutate(
    
    best.ear = ifelse(ear.left < ear.right, yes = ear.left, no = ear.right),
    #ear.diff = ear.left - ear.right, # lower values indicate better hearing
    var = 100 - (50*best.eye),
    var.diff = 50*eye.diff
    
  ) %>%
    
    dplyr::select(
           ID_1,
           id,
           id.twins,
           sex,
           age.weeks.11,
           # hearing (for comparison)
           hear.imp, 
           hear.loss,
           hf.hear.loss,
           best.ear,
           #ear.diff,
           # vision
           var,
           var.diff,
           mother.education,
           asd,
           # cognitive measures
           read.7,
           nonword.rep,
           wold.comp, 
           sum.ccc,
           wisc.total.iq,
           wisc.verbal.iq,
           wisc.performance.iq,
           wisc.arithmetic,
           GCSE
           
           ) %>%
    
  filter(!is.na(var) & !is.na(var.diff))
  
  norm.data$ID_1 <- as.character(norm.data$ID_1)

  return(norm.data)

}

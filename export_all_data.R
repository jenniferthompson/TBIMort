## -- Download TBI data by event -------------------------------------------------------------------
library(RCurl)

## Read in API token from a supersecret file
source('api_info.R')

## Set constant variables
use.uri <- 'https://redcap.vanderbilt.edu/api/'
use.format <- 'csv'
use.rawlabel <- 'label'
use.rawheader <- 'raw'

## -- Export metadata (data dictionary) ------------------------------------------------------------
dict.data <- read.csv(text = postForm(uri = use.uri,
                                      token = my.tbi.token,
                                      content = 'metadata',
                                      format = use.format,
                                      returnFormat = use.format),
                      stringsAsFactors = FALSE)

## -- Read in demographic data - only collected on day 00 ------------------------------------------
demog.data <- read.csv(text = postForm(uri = use.uri,
                                       token = my.tbi.token,
                                       content = 'record',
                                       format = use.format,
                                       type = 'flat',
                                       'forms[0]' = 'demographics_tracs_data',
                                       'events[0]' = 'inpatient_day_00_arm_1',
                                       rawOrLabel = use.rawlabel,
                                       rawOrLabelHeaders = use.rawheader,
                                       returnFormat = use.format),
                       stringsAsFactors = FALSE)

## -- Read in drug infusion data -------------------------------------------------------------------
infusion.data <- read.csv(text = postForm(uri = use.uri,
                                          token = my.tbi.token,
                                          content = 'record',
                                          format = use.format,
                                          type = 'flat',
                                          'fields[0]' = 'mrn',
                                          'forms[0]' = 'drug_infusion_matt_marshall_data',
                                          rawOrLabel = use.rawlabel,
                                          rawOrLabelHeaders = use.rawheader,
                                          returnFormat = use.format),
                          stringsAsFactors = FALSE)

## -- Read in medical record (Ehrenfeld) data ------------------------------------------------------
mr.data <- read.csv(text = postForm(uri = use.uri,
                                    token = my.tbi.token,
                                    content = 'record',
                                    format = use.format,
                                    type = 'flat',
                                    'fields[0]' = 'mrn',
                                    'forms[0]' = 'medical_record_ehrenfeld_data',
                                    rawOrLabel = use.rawlabel,
                                    rawOrLabelHeaders = use.rawheader,
                                    returnFormat = use.format),
                    stringsAsFactors = FALSE)

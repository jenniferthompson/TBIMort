####################################################################################################
## Create analysis data set for TBI/mortality analyses
####################################################################################################

library(Hmisc)
library(dplyr) ## used for easy summarizing for outcome variables

## Function to take a variable with "Missing" as an option and recode as NA
missing.options <- c('Missing', 'Not recorded/Missing', 'Missing/Unknown', 'Unknown/Missing',
                     'Not Recorded', 'Missing/Not Recorded', 'Missing or no data', 'Unknown')
rm.miss <- function(x){ ifelse(x %in% missing.options, NA, x) }

## Function to create and relevel a factor variable given a reference string
f.relevel <- function(x, rlev){ relevel(factor(x), ref = rlev) }

## -- Import data from REDCap ----------------------------------------------------------------------
source('export_all_data.R')

## dict.data: data dictionary
## demog.data: demographics form only (one record per patient)
## mr.data: medical record data (multiple records per patient)
## infusion.data: drug infusion data (multiple records per patient)

## Save typing and change all underscores to dots in variable names
names(demog.data) <- gsub('_', '.', names(demog.data), fixed = TRUE)
names(mr.data) <- gsub('_', '.', names(mr.data), fixed = TRUE)
names(infusion.data) <- gsub('_', '.', names(infusion.data), fixed = TRUE)

## -- Delete rows with completely missing data - artifact of data management error/fix -------------
## These columns can have data: patient ID, REDCap event, "complete" variable
mr.canhave <- c(1, 2, ncol(mr.data))
mr.canthave <- setdiff(1:ncol(mr.data), mr.canhave)
mr.data.miss <- rowSums(!is.na(mr.data[,mr.canthave])) == 0
# mr.data <- mr.data[!mr.data.miss,]

ck.mr.infusion <- merge(mr.data[mr.data.miss, c('mrn', 'redcap.event.name')], infusion.data,
                        by = c('mrn', 'redcap.event.name'),
                        all.x = TRUE, all.y = FALSE)

## -- Merge MR data and drug drip data for calculating drug totals ---------------------------------
daily.data <- merge(mr.data[!mr.data.miss,], infusion.data,
                    by = c('mrn', 'redcap.event.name'), all = TRUE)

## -- Add weight at baseline to daily data for use in SOFA CV, drug calculations -------------------
base.wt.data <- subset(daily.data,
                       redcap.event.name == 'Inpatient Day 00',
                       select = c(mrn, med.wt.encounter))
names(base.wt.data) <- c('mrn', 'base.weight')
daily.data <- merge(daily.data, base.wt.data, by = 'mrn')

## -- Derive necessary variables using daily data --------------------------------------------------
## Pupil reactivity
daily.data$n.pupil.react <-
  ifelse(rowSums(!is.na(daily.data[,c('l.pupil.day', 'r.pupil.day')])) == 0, NA,
         rowSums(daily.data[,c('l.pupil.day', 'r.pupil.day')] == 'Reactive', na.rm = TRUE))
daily.data$pupil.react <- factor(daily.data$n.pupil.react,
                                 levels = 0:2,
                                 labels = c('Both fixed', 'One reactive', 'Both reactive'))

## Mental status:
## coma = minimum RASS = -5 or -4
## delirium = not comatose and CAM+
## normal = not comatose, not CAM+
## missing if RASS missing or > -4 and CAM is missing
daily.data$mental.status <- with(daily.data, {
  factor(ifelse(!is.na(min.rass.score) & min.rass.score %in% c(-4, -5), 3,
         ifelse(is.na(cam.icu), NA,
         ifelse(cam.icu == 'Positive', 2, 1))),
         levels = 1:3, labels = c('Normal', 'Delirious', 'Comatose')) })

## Daily SOFA score ##
## Respiratory component
daily.data$sofa.resp <- with(daily.data, {
  ifelse(is.na(p02.fi02), NA,
  ifelse(p02.fi02 <= 100, 4,
  ifelse(p02.fi02 <= 200, 3,
  ifelse(p02.fi02 <= 300, 2,
  ifelse(p02.fi02 <= 400, 1, 0))))) })

## CNS component
daily.data$sofa.cns <- with(daily.data, {
  ifelse(is.na(min.sum.evm), NA,
  ifelse(min.sum.evm < 6, 4,
  ifelse(min.sum.evm <= 9, 3,
  ifelse(min.sum.evm <= 12, 2,
  ifelse(min.sum.evm <= 14, 1, 0))))) })

## Cardiovascular component
## Create scaled versions of pressors where dosage matters for SOFA: [units]/kg/min
daily.data$scaled.dopa <- with(daily.data, drip.dopa / (base.weight * 24 * 60))
daily.data$scaled.epi <- with(daily.data, drip.epi / (base.weight * 24 * 60))
daily.data$scaled.norepi <- with(daily.data, drip.norepi / (base.weight * 24 * 60))

daily.data$sofa.cv <- with(daily.data, {
  ifelse((!is.na(scaled.dopa) & scaled.dopa > 15) |
           (!is.na(scaled.epi) & scaled.epi > 0.1) |
           (!is.na(scaled.norepi) & scaled.norepi > 0.1), 4,
  ifelse((!is.na(scaled.dopa) & scaled.dopa > 5) |
           (!is.na(drip.epi) & drip.epi > 0) |
           (!is.na(drip.norepi) & drip.norepi > 0), 3,
  ifelse((!is.na(drip.dopa) & drip.dopa > 0) |
           (!is.na(drip.dobu) & drip.dobu > 0) | 
           (!is.na(drip.milri) & drip.milri > 0) |
           (!is.na(drip.vaso) & drip.vaso > 0) |
           (!is.na(drip.phenyleph) & drip.phenyleph > 0), 2,
  ifelse(!is.na(min.map) & min.map < 70, 1,
  ifelse(!is.na(min.map), 0, NA))))) })

## Liver component: per PIs, assume a score of 0 (normal) if no bilirubin available; missingness
##  is highly informative
daily.data$sofa.liver <- with(daily.data, {
  ifelse(is.na(max.tot.bilirubin), 0,
  ifelse(max.tot.bilirubin > 12, 4,
  ifelse(max.tot.bilirubin >= 6, 3,
  ifelse(max.tot.bilirubin >= 2, 2,
  ifelse(max.tot.bilirubin >= 1.2, 1, 0))))) })

## Coagulation component
daily.data$sofa.coag <- with(daily.data, {
  ifelse(is.na(min.platelet), NA,
  ifelse(min.platelet <= 20, 4,
  ifelse(min.platelet <= 50, 3,
  ifelse(min.platelet <= 100, 2,
  ifelse(min.platelet <= 150, 1, 0))))) })

## Renal component
## Double check this once Maddie is done - currently no non-missing values
daily.data$sofa.renal <- with(daily.data, {
  ifelse(is.na(max.creatinine) & is.na(urine), NA,
  ifelse((!is.na(max.creatinine) & max.creatinine > 5) |
           (!is.na(urine) & urine < 200), 4,
  ifelse((!is.na(max.creatinine) & max.creatinine >= 3.5) |
           (!is.na(urine) & urine < 500), 3,
  ifelse(!is.na(max.creatinine) & max.creatinine >= 2, 2,
  ifelse(!is.na(max.creatinine) & max.creatinine >= 1.2, 1, 0))))) })

## Dichotomous version of maximum ICP: normal if score is <= 20, abnormal if > 20
daily.data$max.icp.dich <-
  with(daily.data, factor(ifelse(is.na(max.icp) | max.icp <= 20, 1, 2),
                          levels = 1:2, labels = c('Normal ICP (<=20 or missing)', 'Abnormal (>20)')))

## -- Impute closest value within X days before/after a missing date for several covariates --------
## Create numeric variable for study day
daily.data$study.day <- as.numeric(gsub('Inpatient Day ', '', daily.data$redcap.event.name)) + 1

## Helper function to do imputation for a single-ID data frame
imputer <- function(d1, ## study day or date vector 1 - records that are missing
                    d2, ## study day or date vector 2 - records that are not missing
                    win = 2){ ## number of days to look behind/ahead
  
  if(class(d1) != class(d2)){
    stop('d1 and d2 must be of same class')
  }
  
  ## Create matrix of number of days apart for each missing record vs non-missing records
  ## Columns = each missing row
  ## Rows = distance of each non-missing row from missing row
  if(class(d1) == 'Date'){
    d3 <- sapply(d1, FUN = function(d){ abs(difftime(d, d2, units = 'days')) })
  } else if(class(d1) %in% c('numeric', 'integer')){
    d3 <- sapply(d1, FUN = function(d){abs(d - d2)})
  }
  
  ## If only non-missing day is day 0, sapply returns a vector; force to be a matrix
  if(!is.matrix(d3)){
    d3 <- matrix(d3, nrow = 1)
  }
  
  ## For each column (missing row), get closest non-missing row within "win" days
  apply(d3, 2, function(i) {
    ## Set all non-missing rows > win days away to NA
    i[abs(i) > win] <- NA
    
    if(all(is.na(i))){
      return(NA)
    } else{
      ## Priortize closest, then earliest record
      which.min(i)
    }
  })
}

impute.closest <- function(dataset, ## data set to use - assumes all records are from same ID
                           dayvar,  ## study day/date variable name
                           impvar){ ## name of variable to impute
  
  ## T/F vector indicating whether impvar is missing
  m <- is.na(dataset[, impvar])
  
  ## If any records are missing, replace NA with imputed value
  dataset[,paste0(impvar, '.imp')] <- dataset[,impvar]
  
  if(any(m) & !all(m)){
    ix <- imputer(dataset[m, dayvar], dataset[!m, dayvar], win = 2)
    
    if(any(!is.na(ix))) {
      impvar.imp <- dataset[which(!m)[ix], impvar]
      dataset[m, paste0(impvar, '.imp')] <- dataset[which(!m)[ix], impvar]
    }
  }
  
  dataset
}

## New data frame for each imputed variable
impute.2days <- c('max.motor', 'pupil.react', 'min.glucose', 'min.hemoglobin', 'min.sodium',
                  'sofa.resp', 'sofa.cns', 'sofa.cv', 'sofa.liver', 'sofa.coag', 'sofa.renal')

## List of separate data frames split by ID, sorted by date
daily.data <- daily.data[order(daily.data$mrn, daily.data$study.day),]
daily.data.mrn <- split(daily.data, daily.data$mrn)

## Wrapper function
impute.var <- function(i, impvar){
  tmp <- i[,c('mrn', 'study.day', impvar)]
  impute.closest(dataset = tmp, dayvar = 'study.day', impvar = impvar)
}

## Create a list for each variable imputed using this methodology, containing a list of
##  length(unique(daily.data$mrn)) data frames with mrn, study.day, original and imputed values
imputed.daily <- lapply(impute.2days, FUN = function(impvar){
  lapply(daily.data.mrn, FUN = impute.var, impvar = impvar)
})

imputed.daily.rbind <- lapply(imputed.daily, FUN = function(i){ do.call(rbind, i) })
## rbind is faster
# system.time(imputed.daily.unsplit <- lapply(imputed.daily, FUN = function(i){ unsplit(i, daily.data$mrn) }))

imputed.daily.cbind <- do.call(cbind, lapply(imputed.daily.rbind, FUN = function(i){ i[,ncol(i)] }))
colnames(imputed.daily.cbind) <- paste0(impute.2days, '.imp')

daily.data <- cbind(daily.data, imputed.daily.cbind)

## Function did not preserve levels of pupil.react; re-factor
daily.data$pupil.react.imp <- factor(daily.data$pupil.react.imp - 1,
                                     levels = 0:2,
                                     labels = c('Both fixed', 'One reactive', 'Both reactive'))

## Impute max ICP: if higher than 120, assume 120; if missing, impute uniform value between 0-20
impute.normal <- function(x, min.val, max.val){
  if(!is.na(x)){
    return(x)
  } else{
    return(runif(n = 1, min = min.val, max = max.val))
  }
}

daily.data$max.icp.imp <- apply(as.matrix(daily.data$max.icp),
                                MARGIN = 1,
                                FUN = impute.normal, min.val = 0, max.val = 20)
daily.data$max.icp.imp <- ifelse(daily.data$max.icp.imp > 120, 120, round(daily.data$max.icp.imp))

## Calculate overall and modified SOFA; two versions of each:
## - missing component data assumed to be normal
## - missing component data considered missing
sofa.comps <- paste0('sofa.', c('resp', 'cns', 'cv', 'liver', 'coag', 'renal'))
sofa.mod.comps <- setdiff(sofa.comps, 'sofa.cns')

daily.data$sofa.nanormal <- rowSums(daily.data[,sofa.comps], na.rm = TRUE)
daily.data$sofa.namissing <- rowSums(daily.data[,sofa.comps], na.rm = FALSE)
daily.data$sofa.nanormal.imp <- rowSums(daily.data[,paste0(sofa.comps, '.imp')], na.rm = TRUE)
daily.data$sofa.namissing.imp <- rowSums(daily.data[,paste0(sofa.comps, '.imp')], na.rm = FALSE)

daily.data$sofa.mod.nanormal <- rowSums(daily.data[,sofa.mod.comps], na.rm = TRUE)
daily.data$sofa.mod.namissing <- rowSums(daily.data[,sofa.mod.comps], na.rm = FALSE)
daily.data$sofa.mod.nanormal.imp <- rowSums(daily.data[,paste0(sofa.mod.comps, '.imp')], na.rm = TRUE)
daily.data$sofa.mod.namissing.imp <- rowSums(daily.data[,paste0(sofa.mod.comps, '.imp')], na.rm = FALSE)

## Drug variables - can be done separately from imputation, because we assume no data = dose of 0 ##
## Benzos, in midazolam equivalents
daily.data$bolus.loraz.midaz <- daily.data$bolus.loraz * 2.5
daily.data$drip.loraz.midaz <- daily.data$drip.loraz * 2.5
daily.data$bolus.diaz.midaz <- daily.data$bolus.diaz / 2
daily.data$bolus.alpra.midaz <- daily.data$bolus.alpra / 4
daily.data$bolus.clonaz.midaz <- daily.data$bolus.clonaz / 2

benzo.components <- c('bolus.midaz', 'drip.midaz', 'bolus.loraz.midaz', 'drip.loraz.midaz',
                      'bolus.diaz.midaz', 'bolus.alpra.midaz', 'bolus.clonaz.midaz')

daily.data$tot.benzo <- rowSums(daily.data[,benzo.components], na.rm = TRUE)

## Opioids, in fentanyl equivalents
daily.data$bolus.hydromorph.fent <- (daily.data$bolus.hydromorph / 7.5) * 1000
daily.data$drip.hydromorph.fent <- (daily.data$drip.hydromorph / 7.5) * 1000
daily.data$drip.morph.fent <- (daily.data$drip.morph / 50) * 1000
daily.data$bolus.hydrocod.fent <- (daily.data$bolus.hydrocod * 16.7) / 5
daily.data$bolus.metha.fent <- (daily.data$bolus.metha * 122.1) / 15
daily.data$bolus.remi.fent <- daily.data$bolus.remi / 1.2

opioid.components <- c('bolus.fent', 'drip.fent', 'bolus.hydromorph.fent', 'drip.hydromorph.fent',
                       'drip.morph.fent', 'bolus.hydrocod.fent', 'bolus.metha.fent',
                       'bolus.remi.fent')

daily.data$tot.opioid <- rowSums(daily.data[,opioid.components], na.rm = TRUE)

## Total antipsychotics: use regression formulas in Table 3 of Andreasen et al,
##  Biological Psychiatry, vol 67 issue 3 1 Feb 2010, to get haloperidol equivalents
##  http://www.sciencedirect.com/science/article/pii/S0006322309011251
daily.data$halop.olanz <- (daily.data$bolus.olanz / 2.9)^(1 / 0.805)
daily.data$halop.quet <- (daily.data$bolus.quet / 88.16)^(1 / 0.786)
daily.data$halop.risp <- (daily.data$bolus.risp / 0.79)^(1 / 0.851)
daily.data$halop.zipras <- (daily.data$bolus.zipras / 35.59)^(1 / 0.578)

antipsyc.components <- c('bolus.halop', 'halop.olanz', 'halop.quet', 'halop.risp', 'halop.zipras')

daily.data$tot.antipsyc <- rowSums(daily.data[,antipsyc.components], na.rm = TRUE)

## Total beta blockers, in metoprolol equivalents
daily.data$beta.propran <- daily.data$bolus.propran * 1.25
daily.data$beta.labet <- daily.data$bolus.labet * 0.5
daily.data$beta.esmo <- (daily.data$drip.esmo / 1000) / 10

betablock.components <- c('bolus.metop', 'beta.propran', 'beta.labet', 'beta.esmo')

daily.data$tot.betablock <- rowSums(daily.data[,betablock.components], na.rm = TRUE)

## For single-variable drugs, create new variable that is 0 if no dose noted
daily.data$tot.propofol <- with(daily.data, ifelse(is.na(drip.propofol), 0, drip.propofol))
daily.data$tot.dex <- with(daily.data, ifelse(is.na(drip.dex), 0, drip.dex))
daily.data$tot.pento <- with(daily.data, ifelse(is.na(bolus.pento), 0, bolus.pento))
daily.data$tot.clonid <- with(daily.data, ifelse(is.na(bolus.clonid), 0, bolus.clonid))

## Cube root drug variables
daily.data$tot.opioid.cube <- daily.data$tot.opioid**(1/3)
daily.data$tot.benzo.cube <- daily.data$tot.benzo**(1/3)
daily.data$tot.propofol.cube <- daily.data$tot.propofol**(1/3)
daily.data$tot.dex.cube <- daily.data$tot.dex**(1/3)
daily.data$tot.antipsyc.cube <- daily.data$tot.antipsyc**(1/3)
daily.data$tot.betablock.cube <- daily.data$tot.betablock**(1/3)
daily.data$tot.pento.cube <- daily.data$tot.pento**(1/3)
daily.data$tot.clonid.cube <- daily.data$tot.clonid**(1/3)

## Blood products: convert from mL to units (mL not meaningful clinically); if no data, assume none
daily.data$units.prbc <- with(daily.data, ifelse(is.na(tot.pbrc), 0, tot.pbrc / 350))
daily.data$units.plasma <- with(daily.data, ifelse(is.na(tot.plasma), 0, tot.plasma / 350))
daily.data$units.platelets <- with(daily.data, ifelse(is.na(tot.platelets), 0, tot.platelets / 600))
daily.data$units.cryo <- with(daily.data, ifelse(is.na(tot.cryo), 0, tot.cryo / 250))

## -- Data management for demographic/summary data -------------------------------------------------
demog.data$pt.admit <- as.Date(demog.data$pt.admit, format = '%Y-%m-%d')
demog.data$pt.injury.date <- as.Date(demog.data$pt.injury.date, format = '%Y-%m-%d')
demog.data$hosp.death.date <- as.Date(demog.data$hosp.death.date, format = '%Y-%m-%d')
demog.data$ssdi.death.date <- as.Date(demog.data$ssdi.death.date, format = '%Y-%m-%d')
demog.data$final.death.date <- as.Date(demog.data$final.death.date, format = '%m/%d/%Y')
demog.data$discharge.date <- as.Date(demog.data$discharge.date, format = '%Y-%m-%d')

## -- Calculate outcome variables ------------------------------------------------------------------
## Time to in-hospital death: admission to time of in-hospital death
demog.data$time.death.inhosp <- with(demog.data, {
  ifelse(is.na(pt.admit) | is.na(hosp.death.date), NA,
         as.numeric(difftime(hosp.death.date, pt.admit, units = 'days')) + 1) })

## Time to in-hospital death or discharge: admission to time of in-hospital death, or
##  censored at discharge
demog.data$time.death.dc <- with(demog.data, {
  ifelse(!is.na(time.death.inhosp), time.death.inhosp,
         as.numeric(difftime(discharge.date, pt.admit, units = 'days')) + 1) })

## Time to overall death: admission to max(SSDI death, hospital death)
demog.data$time.death.ever <- with(demog.data, {
  ifelse(is.na(pt.admit) | (is.na(hosp.death.date) & is.na(ssdi.death.date)), NA,
  ifelse(is.na(hosp.death.date) | ssdi.death.date > hosp.death.date,
         as.numeric(difftime(ssdi.death.date, pt.admit, units = 'days')) + 1,
         as.numeric(difftime(hosp.death.date, pt.admit, units = 'days')) + 1)) })

## Time to 3-year death: time to overall death, or censored at 1095 days
demog.data$time.death.3yr <- with(demog.data, {
  ifelse(is.na(time.death.ever) | time.death.ever > 1095, 1095, time.death.ever) })
demog.data$died.3yr <- with(demog.data, {
  factor(ifelse(!is.na(time.death.ever) & time.death.ever <= 1095, 1, 0),
         levels = 0:1, labels = c('No', 'Yes')) })

## Delirium, coma duration
mental.vars <- daily.data %>%
  group_by(mrn) %>%
  summarise(n.recs = n(),
            days.mental = sum(!is.na(mental.status)),
            days.del = ifelse(days.mental == 0, NA,
                              sum(mental.status == 'Delirious', na.rm = TRUE)),
            days.coma = ifelse(days.mental == 0, NA,
                               sum(mental.status == 'Comatose', na.rm = TRUE)))

## DCFDs over the first 14 days after admission
dcfd.data <- daily.data %>%
  filter(redcap.event.name %in%
           paste('Inpatient Day',
                 c('01', '02', '03', '04', '05', '06', '07', '08', '09', '10', '11', '12',
                   '13', '14'))) %>%
  group_by(mrn) %>%
  summarise(days.mental.14 = sum(!is.na(mental.status)),
            dcfd.14 = ifelse(days.mental.14 == 0, NA,
                             14 - sum(mental.status %in% c('Delirious', 'Comatose'),
                                      na.rm = TRUE)))

## -- Demographic data management ------------------------------------------------------------------
demog.data$race2 <- f.relevel(demog.data$race, 'White')
demog.data$gender <- f.relevel(demog.data$gender, 'Male')
demog.data$insurance.code <- f.relevel(demog.data$insurance.code, 'Private')
demog.data$cpr.yn <- f.relevel(demog.data$cpr.yn, 'Not Performed Pre-Hospital, Ref Hospital, or ED')
demog.data$pt.marshall <- f.relevel(demog.data$pt.marshall, 'Marshall Class I')
demog.data$pt.cerebral.na <- f.relevel(rm.miss(demog.data$pt.cerebral), 'No')
demog.data$pt.epidural.na <- f.relevel(rm.miss(demog.data$pt.epidural), 'No')
demog.data$pt.injury.na <- f.relevel(rm.miss(demog.data$pt.injury), 'Blunt')
demog.data$disposition.coded <- f.relevel(demog.data$disposition.coded, 'Discharged home')
demog.data$hosp.death <- f.relevel(demog.data$hosp.death, 'No')
demog.data$final.death <- f.relevel(demog.data$final.death, 'No')

## -- Get baseline MR variables from day 00 --------------------------------------------------------
day00.data <- subset(daily.data,
                     redcap.event.name == 'Inpatient Day 00',
                     select = c(mrn, med.wt.encounter, max.motor, max.motor.imp, min.glucose,
                                min.glucose.imp, min.hemoglobin, min.hemoglobin.imp, pupil.react,
                                pupil.react.imp))
names(day00.data) <- c('mrn', 'base.weight', 'base.motor', 'base.motor.imp', 'base.glucose',
                       'base.glucose.imp', 'base.hemoglobin', 'base.hemoglobin.imp',
                       'base.pupil.react', 'base.pupil.react.imp')

## -- Create analysis data sets --------------------------------------------------------------------
## Use age filter to remove patients with only pupil reactivity data; some sort of weirdness there
tbi.oneobs <- subset(demog.data,
                     select = c(mrn, age, gender, race, insurance.code, iss, cpr.yn, pt.marshall,
                                pt.cerebral, pt.cerebral.na, pt.epidural, pt.epidural.na, pt.injury,
                                pt.injury.na, vent.days, disposition.coded, fim.total, icu.los, los,
                                hosp.death, time.death.inhosp, time.death.dc, died.3yr,
                                time.death.3yr, time.death.ever)) %>%
  left_join(day00.data, by = 'mrn') %>%
  left_join(mental.vars, by = 'mrn') %>%
  left_join(select(dcfd.data, mrn, dcfd.14), by = 'mrn') %>%
  filter(!is.na(age)) %>%
  ## If baseline motor score is still missing after imputation, assume score of 6
  mutate(base.motor.imp = ifelse(is.na(base.motor.imp), 6, base.motor.imp))

label(tbi.oneobs$mrn) <- 'Medical record number'
label(tbi.oneobs$age) <- 'Age at admission'
label(tbi.oneobs$gender) <- 'Gender'
label(tbi.oneobs$race) <- 'Race'
label(tbi.oneobs$insurance.code) <- 'Insurance type'
label(tbi.oneobs$iss) <- 'Injury Severity Score'
label(tbi.oneobs$cpr.yn) <- 'CPR performed'
label(tbi.oneobs$pt.marshall) <- 'Marshall Class Equivalent'
label(tbi.oneobs$pt.cerebral) <- 'Cerebral subarachnoid hemorrhage'
label(tbi.oneobs$pt.cerebral.na) <- 'Cerebral SAH (unknown = NA)'
label(tbi.oneobs$pt.epidural) <- 'Cerebral extra/epidural mass'
label(tbi.oneobs$pt.epidural.na) <- 'Cerebral extra/epidural mass (unknown = NA)'
label(tbi.oneobs$pt.injury) <- 'Injury type'
label(tbi.oneobs$pt.injury.na) <- 'Injury type (unknown = NA)'
label(tbi.oneobs$vent.days) <- 'Ventilator days'
label(tbi.oneobs$disposition.coded) <- 'Discharge disposition'
label(tbi.oneobs$fim.total) <- 'FIM total score'
label(tbi.oneobs$icu.los) <- 'ICU LOS'
label(tbi.oneobs$los) <- 'Hospital LOS'
label(tbi.oneobs$hosp.death) <- 'Died in hospital'
label(tbi.oneobs$time.death.inhosp) <- 'Days to in-hospital death'
label(tbi.oneobs$time.death.dc) <- 'Hospital LOS (days to in-hospital death or discharge)'
label(tbi.oneobs$time.death.ever) <- 'Days to death'
label(tbi.oneobs$died.3yr) <- 'Died on or before 3-year mark'
label(tbi.oneobs$time.death.3yr) <- 'Days to death or 3-year mark'
label(tbi.oneobs$base.weight) <- 'Median weight during encounter (kg)'
label(tbi.oneobs$base.motor) <- 'Maximum motor response, day 0'
label(tbi.oneobs$base.motor.imp) <- 'Max motor response, day 0 (imputed)'
label(tbi.oneobs$base.glucose) <- 'Minimum glucose, day 0'
label(tbi.oneobs$base.glucose.imp) <- 'Min glucose, day 0 (imputed)'
label(tbi.oneobs$base.hemoglobin) <- 'Minimum hemoglobin, day 0'
label(tbi.oneobs$base.hemoglobin.imp) <- 'Min hemoglobin, day 0 (imputed)'
label(tbi.oneobs$base.pupil.react) <- 'Pupil reactivity, day 0'
label(tbi.oneobs$base.pupil.react.imp) <- 'Pupil reactivity, day 0 (imputed)'
label(tbi.oneobs$n.recs) <- 'Number of daily records'
label(tbi.oneobs$days.mental) <- 'Days with mental status info'
label(tbi.oneobs$days.del) <- 'Days of delirium'
label(tbi.oneobs$days.coma) <- 'Days of coma'
label(tbi.oneobs$dcfd.14) <- 'DCFDs within 14 days of admission'

tbi.daily <- subset(daily.data,
                    mrn %in% tbi.oneobs$mrn,
                    select = c(mrn, redcap.event.name, study.day, mental.status, max.motor,
                               max.motor.imp, pupil.react, pupil.react.imp, min.glucose,
                               min.glucose.imp, min.hemoglobin, min.hemoglobin.imp, min.sodium,
                               min.sodium.imp, max.icp, max.icp.imp, max.icp.dich,
                               sofa.resp, sofa.resp.imp, sofa.cns, sofa.cns.imp, sofa.cv,
                               sofa.cv.imp, sofa.liver, sofa.liver.imp, sofa.coag, sofa.coag.imp,
                               sofa.renal, sofa.renal.imp, sofa.nanormal, sofa.nanormal.imp,
                               sofa.namissing, sofa.namissing.imp, sofa.mod.nanormal,
                               sofa.mod.nanormal.imp, sofa.mod.namissing, sofa.mod.namissing.imp,
                               tot.benzo, tot.benzo.cube, tot.opioid, tot.opioid.cube, tot.propofol,
                               tot.propofol.cube, tot.dex, tot.dex.cube, tot.antipsyc,
                               tot.antipsyc.cube, tot.betablock, tot.betablock.cube, tot.pento,
                               tot.pento.cube, tot.clonid, tot.clonid.cube, units.cryo, units.plasma,
                               units.platelets, units.prbc))

names(tbi.daily) <- gsub('^redcap\\.event\\.name$', 'event', names(tbi.daily))

## Remove records from tbi.daily which are after in-hospital death or discharge date
tbi.daily <- tbi.daily %>%
  left_join(dplyr::select(tbi.oneobs, mrn, time.death.dc), by = 'mrn') %>%
  filter(study.day <= time.death.dc) %>%
  dplyr::select(-time.death.dc)

label(tbi.daily$mrn) <- 'Medical record number'
label(tbi.daily$event) <- 'Study event'
label(tbi.daily$study.day) <- 'Study day'
label(tbi.daily$mental.status) <- 'Mental status'
label(tbi.daily$max.motor) <- 'Maximum motor response'
label(tbi.daily$max.motor.imp) <- 'Max motor response (imputed)'
label(tbi.daily$pupil.react) <- 'Pupil reactivity'
label(tbi.daily$pupil.react.imp) <- 'Pupil reactivity (imputed)'
label(tbi.daily$min.glucose) <- 'Minimum glucose'
label(tbi.daily$min.glucose.imp) <- 'Min glucose (imputed)'
label(tbi.daily$min.hemoglobin) <- 'Minimum hemoglobin'
label(tbi.daily$min.hemoglobin.imp) <- 'Min hemoglobin (imputed)'
label(tbi.daily$min.sodium) <- 'Minimum sodium'
label(tbi.daily$min.sodium.imp) <- 'Min sodium (imputed)'
label(tbi.daily$max.icp) <- 'Maximum ICP'
label(tbi.daily$max.icp.imp) <- 'Max ICP (imputed)'
label(tbi.daily$max.icp.dich) <- 'Max ICP, dichotomous'
label(tbi.daily$sofa.resp) <- 'Respiratory SOFA'
label(tbi.daily$sofa.cns) <- 'CNS SOFA'
label(tbi.daily$sofa.cv) <- 'Cardiovascular SOFA'
label(tbi.daily$sofa.liver) <- 'Liver SOFA'
label(tbi.daily$sofa.coag) <- 'Coagulation SOFA'
label(tbi.daily$sofa.renal) <- 'Renal SOFA'
label(tbi.daily$sofa.resp.imp) <- 'Respiratory SOFA (imputed)'
label(tbi.daily$sofa.cns.imp) <- 'CNS SOFA (imputed)'
label(tbi.daily$sofa.cv.imp) <- 'Cardiovascular SOFA (imputed)'
label(tbi.daily$sofa.liver.imp) <- 'Liver SOFA (imputed)'
label(tbi.daily$sofa.coag.imp) <- 'Coagulation SOFA (imputed)'
label(tbi.daily$sofa.renal.imp) <- 'Renal SOFA (imputed)'
label(tbi.daily$sofa.nanormal) <- 'Overall SOFA, missing values considered normal'
label(tbi.daily$sofa.namissing) <- 'Overall SOFA, missing values left as missing'
label(tbi.daily$sofa.mod.nanormal) <- 'Modified SOFA, missing values considered normal'
label(tbi.daily$sofa.mod.namissing) <- 'Modified SOFA, missing values left as missing'
label(tbi.daily$tot.benzo) <- '24h benzodiazepines, midaz. equivalents'
label(tbi.daily$tot.benzo.cube) <- '24h benzodiazepines, cube root'
label(tbi.daily$tot.opioid) <- '24h opioids, fentanyl equivalents'
label(tbi.daily$tot.opioid.cube) <- '24h opioids, cube root'
label(tbi.daily$tot.propofol) <- '24h propofol'
label(tbi.daily$tot.propofol.cube) <- '24h propofol, cube root'
label(tbi.daily$tot.dex) <- '24h dexmedetomidine'
label(tbi.daily$tot.dex.cube) <- '24h dexmedetomidine, cube root'
label(tbi.daily$tot.antipsyc) <- '24h antipsychotics, haloperidol equivalents'
label(tbi.daily$tot.antipsyc.cube) <- '24h antipsychotics, cube root'
label(tbi.daily$tot.betablock) <- '24h beta blockers'
label(tbi.daily$tot.betablock.cube) <- '24h beta blockers, cube root'
label(tbi.daily$tot.pento) <- '24h pentobarbital'
label(tbi.daily$tot.pento.cube) <- '24h pentobarbital, cube root'
label(tbi.daily$tot.clonid) <- '24h clonidine'
label(tbi.daily$tot.clonid.cube) <- '24h clonidine, cube root'
label(tbi.daily$units.prbc) <- 'Units of packed red blood cells given (mL/350)'
label(tbi.daily$units.plasma) <- 'Units of plasma given (mL/350)'
label(tbi.daily$units.platelets) <- 'Units of platelets given (mL/600)'
label(tbi.daily$units.cryo) <- 'Units of cryoprecipitate given (mL/250)'

## Save date that analysis data sets were created
datasets.created.at <- Sys.time()

save(tbi.oneobs, tbi.daily, datasets.created.at, file = 'AnalysisData/tbi_datasets.Rdata')

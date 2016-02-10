####################################################################################################
## Create analysis data set for TBI/mortality analyses
####################################################################################################

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

## -- Merge MR data and drug drip data for calculating drug totals ---------------------------------
daily.data <- merge(mr.data, infusion.data, by = c('mrn', 'redcap.event.name'), all = TRUE)

## -- Derive necessary variables using daily data --------------------------------------------------
## Pupil reactivity
## Double check this once Maddie is done entering data; currently left pupil has 23% missingness
##  but right pupil has nearly 99%!
daily.data$n.pupil.react <-
  ifelse(rowSums(!is.na(daily.data[,c('l.pupil.day', 'r.pupil.day')])) == 0, NA,
         rowSums(daily.data[,c('l.pupil.day', 'r.pupil.day')] == 'Reactive', na.rm = TRUE))
daily.data$pupil.react <- factor(daily.data$n.pupil.react,
                                 levels = 0:2,
                                 labels = c('Both fixed', 'One reactive', 'Both reactive'))

## -- Get baseline MR variables from day 00 (need pupil reactivity) --------------------------------
day00.data <- subset(daily.data,
                     redcap.event.name == 'Inpatient Day 00',
                     select = c(mrn, med.wt.encounter, max.motor, min.glucose, min.hemoglobin,
                                pupil.react))
names(day00.data) <- c('mrn', 'base.weight', 'base.motor', 'base.glucose', 'base.hemoglobin',
                       'base.pupil.react')

## Merge baseline weight back onto daily data for use in drug scaling
daily.data <- merge(daily.data,
                    subset(day00.data, select = c(mrn, base.weight)),
                    by = 'mrn', all.x = TRUE, all.y = FALSE)

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
## Double check this once Maddie is done - currently no non-missing values
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

## Calculate overall and modified SOFA; two versions of each:
## - missing component data assumed to be normal
## - missing component data considered missing
sofa.comps <- paste0('sofa.', c('resp', 'cns', 'cv', 'liver', 'coag', 'renal'))
sofa.mod.comps <- setdiff(sofa.comps, 'sofa.cns')

daily.data$sofa.nanormal <- rowSums(daily.data[,sofa.comps], na.rm = TRUE)
daily.data$sofa.namissing <- rowSums(daily.data[,sofa.comps], na.rm = FALSE)

daily.data$sofa.mod.nanormal <- rowSums(daily.data[,sofa.mod.comps], na.rm = TRUE)
daily.data$sofa.mod.namissing <- rowSums(daily.data[,sofa.mod.comps], na.rm = FALSE)

## -- Data management for demographic/summary data -------------------------------------------------
demog.data$pt.admit <- as.Date(demog.data$pt.admit, format = '%Y-%m-%d')
demog.data$pt.injury.date <- as.Date(demog.data$pt.injury.date, format = '%Y-%m-%d')
demog.data$hosp.death.date <- as.Date(demog.data$hosp.death.date, format = '%Y-%m-%d')
demog.data$ssdi.death.date <- as.Date(demog.data$ssdi.death.date, format = '%Y-%m-%d')
demog.data$final.death.date <- as.Date(demog.data$final.death.date, format = '%m/%d/%Y')
demog.data$discharge.date <- as.Date(demog.data$discharge.date, format = '%Y-%m-%d')

## -- Calculate outcome variables ------------------------------------------------------------------
## Time to in-hospital death: admission to time of in-hospital death, or censored at discharge
demog.data$time.death.inhosp <- with(demog.data, {
  ifelse(is.na(pt.admit), NA,
  ifelse(!is.na(hosp.death.date),
         as.numeric(difftime(hosp.death.date, pt.admit, units = 'days')),
         as.numeric(difftime(discharge.date, pt.admit, units = 'days')))) })

## Time to overall death: admission to max(SSDI death, hospital death)
demog.data$time.death.ever <- with(demog.data, {
  ifelse(is.na(pt.admit) | (is.na(hosp.death.date) & is.na(ssdi.death.date)), NA,
  ifelse(is.na(hosp.death.date) | ssdi.death.date > hosp.death.date,
         as.numeric(difftime(ssdi.death.date, pt.admit, units = 'days')),
         as.numeric(difftime(hosp.death.date, pt.admit, units = 'days')))) })

## Time to 3-year death: time to overall death, or censored at 1095 days
demog.data$time.death.3yr <-
  ifelse(is.na(demog.data$time.death.ever), 1095, demog.data$time.death.ever)

## Delirium, coma duration
library(dplyr)
mental.vars <- daily.data %>%
  group_by(mrn) %>%
  summarise(n.recs = n(),
            days.assessed = sum(!is.na(mental.status)),
            days.del = ifelse(days.assessed == 0, NA,
                              sum(mental.status == 'Delirious', na.rm = TRUE)),
            days.coma = ifelse(days.assessed == 0, NA,
                               sum(mental.status == 'Comatose', na.rm = TRUE)))

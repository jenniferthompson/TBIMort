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

## -- Derive necessary variables -------------------------------------------------------------------
## Pupil reactivity
## Double check this once Maddie is done entering data; currently left pupil has 23% missingness
##  but right pupil has nearly 99%!
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
## Have question out to M/J about pressor doses
daily.data$sofa.cv <- with(daily.data, {
  ifelse((!is.na(drip.dopa) & drip.dopa > 15) |
           (!is.na(drip.epi) & drip.epi > 0.1) |
           (!is.na(drip.norepi) & drip.norepi > 0.1), 4,
  ifelse((!is.na(drip.dopa) & drip.dopa > 5) |
           (!is.na(drip.epi) & drip.epi > 0) |
           (!is.na(drip.norepi) & drip.norepi > 0), 3,
  ifelse((!is.na(drip.dopa) & drip.dopa > 0) |
           (!is.na(drip.dobu) & drip.dobu > 0) | 
           (!is.na(drip.milri) & drip.milri > 0) |
           (!is.na(drip.vaso) & drip.vaso > 0) |
           (!is.na(drip.phenyleph) & drip.phenyleph > 0), 2,
  ifelse(!is.na(min.map) & min.map < 70, 1,
  ifelse(!is.na(min.map), 0, NA))))) })

## Liver component
## Double check this once Maddie is done - currently no non-missing values
daily.data$sofa.liver <- with(daily.data, {
  ifelse(is.na(max.tot.bilirubin), NA,
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

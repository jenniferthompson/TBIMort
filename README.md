# TBIMort

This project uses EMR and SSDI data to look at risk factors for long-term mortality and other
secondary outcomes among traumatic brain injury patients.

R Scripts
---------
* `api_info.R`: Contains API token for exporting data.
* [`export_all_data.R`](https://github.com/jenniferthompson/TBIMort/blob/master/export_all_data.R): Exports raw data from a REDCap database.
* [`tbimortality_datamgmt.R`](https://github.com/jenniferthompson/TBIMort/blob/master/tbimortality_datamgmt.R): Creates analysis data sets `tbi.oneobs` and `tbi.daily` from raw data and saves in /AnalysisData.
* [`tbi_timevarying.R`](https://github.com/jenniferthompson/TBIMort/blob/master/tbi_timevarying.R): Function to create data sets for time-varying Cox models.

Rnw Scripts
-----------
* [`tbi_descstats.Rnw`](https://github.com/jenniferthompson/TBIMort/blob/master/tbi_descstats.Rnw): Produces descriptive statistics for variables in `tbi.oneobs` and `tbi.daily`, both detailed (using the `Hmisc` `describe()` function) and summarized, including missingness both for single records and for patients for daily variables.
* [`tbi_mental_nomental.Rnw`](https://github.com/jenniferthompson/TBIMort/blob/master/tbi_mental_nomental.Rnw): Describes patients who had at least one CAM/RASS assessment recorded vs those who did not.
* [`tbi_sccm2017_abstract.Rnw`](https://github.com/jenniferthompson/TBIMort/blob/master/tbi_sccm2017_abstract.Rnw): Preliminary analyses for submission to SCCM 2017 conference.

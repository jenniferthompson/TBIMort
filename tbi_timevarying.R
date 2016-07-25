## -- Function to create time-varying Cox data -- ##
## -- Original data set requires "id", "study.day" variables -- ##
timevarying.data <-
  function(idvar = 'id', ## string; name of cluster ID variable
           recvar = 'study.day', ## string; name of record indicator; variable must be numeric
           org.data,     ## Original data with multiple records per patient
           time.var,     ## Variable indicating end of records (eg, day of death)
           event.var,    ## Variable indicating whether event happened
           event.string, ## Level of event variable indicating whether event happened
           data.set,     ## Data set (one record/pt) containing time/event variables
           out.strings = c('Alive through end of interval',
                           'Died at end of interval')){ ## Labels for 0/1 values of outcome variable
    
    ## tbl_df objects don't work here; force org.data and data.sets to be data.frames
    org.data <- as.data.frame(org.data)
    data.set <- as.data.frame(data.set)
    
    ## Get vector of unique IDs, create null list with one slot for each
    ids <- unique(org.data[, idvar])
    final.data <- vector("list", length(ids))
    cols2drop <- match(c(idvar, recvar), names(org.data))
    
    ## For each ID...
    for(i in seq_along(ids)) {
      cur.id <- ids[i]
      
      ## Get all relevant daily data for specific patient
      use.data <- org.data[org.data[,idvar] == cur.id, -cols2drop]
      
      ## find rows with missing data
      valid.rows <- apply(use.data, MARGIN = 1, FUN=function(i) !all(is.na(i)))
      
      ## If ID has at least one row with non-missing data, create time-varying data set;
      ## otherwise, return NULL
      if(sum(valid.rows) > 0){
        ## Create vector of study days with non-missing data and first day with non-missing data
        study.day.id <- org.data[org.data[,idvar] == cur.id, recvar][valid.rows]
        first.day <- min(study.day.id)
        
        ## Delete any rows with missing data
        use.data <- use.data[valid.rows,]
        
        ## For each ID, cycle through each record. If current row = last row, increment counter by one.
        nr <- nrow(use.data)
        if(nr > 1) {
          stop.day <- rep(NA, nrow(use.data))
          tmp <- apply(use.data, MARGIN = 1, paste, collapse = '')
          for(j in seq(nr - 1)) {
            if(tmp[j] != tmp[j+1]) {
              stop.day[j] <- study.day.id[j]
            }
          }
        } else {
          stop.day <- NA
        }
        
        ## Cbind each stop day to daily data set, and add date of death/end of period to last row
        stop.day[length(stop.day)] <- data.set[data.set[,idvar] == cur.id, time.var]
        use.data <- cbind(use.data, stop.day)
        
        ## Subset to only non-repeated rows, create start date vector
        use.data <- use.data[!is.na(stop.day),]
        use.data$start.day <- c((first.day - 1), use.data$stop.day[-nrow(use.data)])
        nr <- nrow(use.data)
        
        ## If patient died or discharged on a day with assessments, or on the day of enrollment,
        ##  assign interval half a day
        use.data$stop.day <- ifelse(use.data$stop.day == use.data$start.day,
                                    use.data$stop.day + 0.5,
                                    use.data$stop.day)
        
        ## Create event vector: all 0 unless patient experienced event, in which case last value = 1
        use.data$had.event <- 0
        use.data$had.event[nr] <- (data.set[data.set[,idvar] == cur.id, event.var] == event.string)*1
        
        final.data[[i]] <- cbind(cur.id, use.data)
        
        cat(sprintf("\rFinished %s", cur.id))
      } else{
        final.data[[i]] <- NULL
      }
    }
    cat("\n")
    
    ## combine list of data.sets into one
    final.data <- do.call('rbind', final.data)
    
    final.data$had.event <- factor(final.data$had.event, levels = 0:1, labels = out.strings)
    names(final.data) <- gsub('cur.id', idvar, names(final.data), fixed = TRUE)
    return(final.data)
  }

basemodelcols <- names(dat)[c(1, 10, 14:16, 18:23)]

get.all.scc <- function(dat) {
    basecodes <- sapply(1:nrow(dat), function(ii) paste(dat[ii, basemodelcols], collapse=', '))

    basedat <- data.frame(basecode=basecodes,
                          scc=dat$`Reported Base Model SCC (if applicable)`)

    for (col in names(dat)) {
        if (col == "Central Value ($ per ton CO2)")
            next # Already recorded
        else if (col == "PAPER LOCATION")
            basedat[, col] <- "Reported Base Model SCC" # Clarified
        else if (col %in% names(dat)[(which(names(dat) == "NOTES")+1):length(names(dat))])
            basedat[, col] <- dat[, col] # Added columns: keep values
        else if (col %in% c("Bibtex Name", "Reference"))
            basedat[, col] <- paste(dat[, col], "-base") # Paper-specific reference
        else if (col %in% c("Added By", basemodelcols))
            basedat[, col] <- dat[, col] # Keep values
        else
            basedat[, col] <- NA # Drop values
    }
    basedat$modified <- F

    basedat <- basedat[!duplicated(basedat),]

    df <- rbind(basedat,
                cbind(basecode=basecodes, scc=dat$`Central Value ($ per ton CO2)`,
                      dat[, names(dat) != 'Central Value ($ per ton CO2)'], modified=T))

    for (rr in which(df$`Tipping Points2` == "-1.0")) {
        df$`Tipping Points2`[df$basecode == df$basecode[rr] & !df$modified] <- "1.0"
        df$`Tipping Points2`[rr] <- NA
    }

    df
}

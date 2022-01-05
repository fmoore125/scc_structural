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

multivar.prep <- function(df) {
    names(df)[grep("Alternative ethical approaches", names(df))] <- "Alternative ethical approaches"

    df$`IAM Calibrated To (if applicable)`[is.na(df$`IAM Calibrated To (if applicable)`)] <- "None"
    df$`Backstop Price?`[is.na(df$`Backstop Price?`)] <- "0"
    df$`Other Market Failure?`[is.na(df$`Other Market Failure?`)] <- "0"
    df$`Other Market Failure?`[df$`Other Market Failure?` != "0"] <- "1.0"
    df$`Market Only Damages`[is.na(df$`Market Only Damages`)] <- "0"
    df$`Carbon Cycle`[is.na(df$`Carbon Cycle`)] <- "0"
    df$`Climate Model`[is.na(df$`Climate Model`)] <- "0"
    df$`Tipping Points`[is.na(df$`Tipping Points`)] <- "0"
    df$`Tipping Points2`[is.na(df$`Tipping Points2`)] <- "0"
    df$`Persistent / Growth Damages`[is.na(df$`Persistent / Growth Damages`)] <- "0"
    df$`Epstein-Zin`[is.na(df$`Epstein-Zin`)] <- "0"
    df$`Ambiguity/Model Uncertainty`[is.na(df$`Ambiguity/Model Uncertainty`)] <- "0"
    df$`Limitedly-Substitutable Goods`[is.na(df$`Limitedly-Substitutable Goods`)] <- "0"
    df$`Inequality Aversion`[is.na(df$`Inequality Aversion`)] <- "0"
    df$`Learning`[is.na(df$`Learning`)] <- "0"
    df$`Alternative ethical approaches`[is.na(df$`Alternative ethical approaches`)] <- "0"

    df
}

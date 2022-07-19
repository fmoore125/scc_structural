## setwd("~/research/scciams/scc_structural")

library(readxl)
library(cluster)

df <- read_excel("../data/data_SCC-expert-survey_final_anonymous.xlsx")
df$anon <- df$Anonymous
df$anon[is.na(df$anon)] <- 1 - df$Non_anonymous[is.na(df$anon)]

used <- df[!is.na(df$SCC_Lit_Central) & !is.na(df$SCC_True_Central),
           c('SCC_Lit_2p5perc', 'SCC_Lit_Central', 'SCC_Lit97p5perc', 'SCC_True_2p5perc', 'SCC_True_Central', 'SCC_True_97p5perc',
             'Q2_MeanSCC_earth_system', 'Q2_MeanSCC_tipping_climate_system', 'Q2_MeanSCC_tipping_damages', 'Q2_MeanSCC_persistent_damages',
             'Q2_MeanSCC_substitutability', 'Q2_MeanSCC_Epstein-Zin', 'Q2_MeanSCC_model_uncertainty', 'Q2_MeanSCC_learning',
             'Q2_MeanSCC_distributional_weighting', 'Q2_Include_earth_system', 'Q2_Include_tipping_climate_system',
             'Q2_Include_tipping_damages', 'Q2_Include_persistent_damage', 'Q2_Include_substitutability', 'Q2_Include_Epstein-Zin',
             'Q2_Include_model_uncertainty', 'Q2_Include_learning', 'Q2_Include_distributional_weighting', 'SCC_wedge_value',
             'SCC_wedge', 'SCCwedge_Drivers_earth_system', 'SCCwedge_Drivers_tipping_damage', 'SCCwedge_Drivers_persistent_damages',
             'SCCwedge_Drivers_substitutability', 'SCCwedge_Drivers_Epstein-Zin', 'SCCwedge_Drivers_model_uncertainty',
             'SCCwedge_Drivers_learning', 'SCCwedge_Drivers_distributional_weighting', 'SCCwedge_Drivers_endogenous_tech_progress',
             'SCCwedge_Drivers_ethical_approaches', 'SCCwedge_Drivers_pure_time_pref', 'SCCwedge_Drivers_EMUC',
             'SCCwedge_Drivers_damage_function_parameters', 'SCCwedge_Drivers_adaptation', 'SCCwedge_Drivers_tipping_climate_system',
             'SCCwedge_Drivers_Other', 'anon')]
mat <- used
for (col in names(mat)) {
    values <- mat[, col, drop=T]
    mat[, col] <- (values - mean(values, na.rm=T)) / sd(values, na.rm=T)
}

distmat <- matrix(NA, nrow(mat), nrow(mat))
colnames(distmat) <- 1:nrow(mat)
rownames(distmat) <- 1:nrow(mat)

for (ii in 1:(nrow(mat)-1)) {
    for (jj in (ii+1):nrow(mat)) {
        vali <- t(mat[ii,])
        valj <- t(mat[jj,])
        both <- !is.na(vali) & !is.na(valj)
        none <- is.na(vali) & is.na(valj)
        distmat[jj, ii] <- sum((vali[both] - valj[both])^2) + sum(!both) - sum(none)
    }
}

clusters <- hclust(as.dist(distmat))
plot(clusters)

ccut <- cutree(clusters, 6)

clcore <- colMeans(used[ccut == 5,], na.rm=T)
clleft <- colMeans(used[ccut == 1,], na.rm=T)
clright <- colMeans(used[!(ccut %in% c(1, 5)),], na.rm=T)

sum(clcore[27:42])
sum(clleft[27:42])
sum(clright[27:42])

rbind(c(count=sum(ccut == 5), clcore),
      c(count=sum(ccut == 1), clleft),
      c(count=sum(!(ccut %in% c(1, 5))), clright))


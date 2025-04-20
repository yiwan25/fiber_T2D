# ceraf, frtaf, vegaf, stvegaf, nonstvegaf, all mutually adjusted

library(tidyverse)
library(dplyr)

nhs1 <- read.csv("/udd/nhyiw/Fiber/biomarker/readnhs_for_bld.csv", header = T)
nhs2 <- read.csv("/udd/nhyiw/Fiber/biomarker/readnhs2_for_bld.csv", header = T)
hpfs <- read.csv("/udd/nhyiw/Fiber/biomarker/readhpfs_for_bld.csv", header = T)

nhs1$cohort="nh"
nhs2$cohort="n2"
hpfs$cohort="hp"

pool=bind_rows(nhs1, nhs2)
pool=bind_rows(pool, hpfs)
pool <- pool %>% filter(canbase==0 & hrtbase==0 & strbase==0 & dbbase==0 & is.na(dbbase2))

##############################################
#       fiber-biomarker associations
##############################################
pool$hp <- ifelse(pool$cohort=="hp",1,0)
pool$n2 <- ifelse(pool$cohort=="n2",1,0)

pool[pool$cohort=="hp",]$pmh <- 9

bio <- c("adj_lncpep", "adj_lna1c", "adj_lnchol",  "adj_lnhdl",  "adj_lnldl",  "adj_lntrig", "adj_lnratio",
         "adj_lncrp", "adj_lnadipo", "adj_lnleptin")

results_aofib <- data.frame(type="aofib",
                           bio=bio,
                           n=NA,
                           beta=NA,
                           lci=NA,
                           uci=NA,
                           p=NA)
results_ceraf <- data.frame(type="ceraf",
                              bio=bio,
                              n=NA,
                              beta=NA,
                              lci=NA,
                              uci=NA,
                              p=NA)
results_vegaf <- data.frame(type="vegaf",
                           bio=bio,
                           n=NA,
                           beta=NA,
                           lci=NA,
                           uci=NA,
                           p=NA)
results_frtaf <- data.frame(type="frtaf",
                              bio=bio,
                              n=NA,
                              beta=NA,
                              lci=NA,
                              uci=NA,
                              p=NA)
results_stvegaf <- data.frame(type="stvegaf",
                              bio=bio,
                              n=NA,
                              beta=NA,
                              lci=NA,
                              uci=NA,
                              p=NA)
results_nonstvegaf <- data.frame(type="nonstvegaf",
                              bio=bio,
                              n=NA,
                              beta=NA,
                              lci=NA,
                              uci=NA,
                              p=NA)

pool <- pool %>% mutate(bmi_avg=ifelse(is.na(bmi_avg), median(bmi_avg, na.rm=T), bmi_avg),
                        pa_avg=ifelse(is.na(pa_avg), median(pa_avg, na.rm=T), pa_avg),
                        trans_avg=ifelse(is.na(trans_avg), median(trans_avg, na.rm=T), trans_avg),
                        nses_avg=ifelse(is.na(nses_avg), median(nses_avg, na.rm=T), nses_avg),
                        asp=ifelse(asp==9|is.na(asp), 0, asp),
                        pmh=ifelse(is.na(pmh), 1, pmh),
                        cigg=case_when(is.na(cigg) & pkyrg==1 ~ 1,
                                       is.na(cigg) & pkyrg!=1 ~ 2,
                                       TRUE ~ cigg))

check_na <- data.frame(var=colnames(pool),
                       na=colSums(is.na(pool)))

pool <- pool %>% mutate(pmh=factor(pmh),
                        asp=factor(asp),
                        cigg=factor(cigg),
                        pkyrg=factor(pkyrg))

quantile(pool$aofib_avg, probs = 0.9) - quantile(pool$aofib_avg, probs = 0.1) # 13.345
quantile(pool$ceraf_avg, probs = 0.9) - quantile(pool$ceraf_avg, probs = 0.1) # 6.52
quantile(pool$frtaf_avg, probs = 0.9) - quantile(pool$frtaf_avg, probs = 0.1) # 5.555
quantile(pool$vegaf_avg, probs = 0.9) - quantile(pool$vegaf_avg, probs = 0.1) # 6.5
quantile(pool$stvegaf_avg, probs = 0.9) - quantile(pool$stvegaf_avg, probs = 0.1) # 1.475694
quantile(pool$nonstvegaf_avg, probs = 0.9) - quantile(pool$nonstvegaf_avg, probs = 0.1) # 5.85

pool$aofib_avg <- pool$aofib_avg/10
pool$ceraf_avg <- pool$ceraf_avg/5
pool$frtaf_avg <- pool$frtaf_avg/5
pool$vegaf_avg <- pool$vegaf_avg/5
pool$nonstvegaf_avg <- pool$nonstvegaf_avg/5
pool$stvegaf_avg <- pool$stvegaf_avg/2

for (i in 1:10) {
  tempdata <- pool %>% filter(!is.na(aofib_avg))
  tempdata <- tempdata[!is.na(tempdata[,bio[i]]),]

  # aofib
  tempmod1 <- lm(paste0(bio[i],"~aofib_avg+age_blddraw+race+fast+caco+n2+pmh+asp+mv+hbpbase+cholbase+nses_avg+
  bmi_avg+pa_avg+cigg+pkyrg+alco_avg+calor_avg+poly_avg+sat_avg+trans_avg+gl_avg"), data=tempdata)
  tempsum1 <- summary(tempmod1)
  results_aofib[i, "n"] <- dim(tempdata)[1]
  results_aofib[i, "beta"] <- (exp(tempsum1$coefficients[2,1])-1)*100
  results_aofib[i, "lci"] <- (exp(confint(tempmod1)[2,1])-1)*100
  results_aofib[i, "uci"] <- (exp(confint(tempmod1)[2,2])-1)*100
  results_aofib[i, "p"] <- tempsum1$coefficients[2,4]
  
  # ceraf
  tempmod2 <- lm(paste0(bio[i],"~ceraf_avg+frtaf_avg+vegaf_avg+age_blddraw+race+fast+caco+n2+pmh+asp+mv+nses_avg+
  hbpbase+cholbase+bmi_avg+pa_avg+cigg+pkyrg+alco_avg+calor_avg+poly_avg+sat_avg+trans_avg+gl_avg"), data=tempdata)
  tempsum2 <- summary(tempmod2)
  results_ceraf[i, "n"] <- dim(tempdata)[1]
  results_ceraf[i, "beta"] <- (exp(tempsum2$coefficients[2,1])-1)*100
  results_ceraf[i, "lci"] <- (exp(confint(tempmod2)[2,1])-1)*100
  results_ceraf[i, "uci"] <- (exp(confint(tempmod2)[2,2])-1)*100
  results_ceraf[i, "p"] <- tempsum2$coefficients[2,4]
  
  # vegaf
  tempmod3 <- lm(paste0(bio[i],"~vegaf_avg+ceraf_avg+frtaf_avg+age_blddraw+race+fast+caco+n2+pmh+asp+mv+nses_avg+
  hbpbase+cholbase+bmi_avg+pa_avg+cigg+pkyrg+alco_avg+calor_avg+poly_avg+sat_avg+trans_avg+gl_avg"), data=tempdata)
  tempsum3 <- summary(tempmod3)
  results_vegaf[i, "n"] <- dim(tempdata)[1]
  results_vegaf[i, "beta"] <- (exp(tempsum3$coefficients[2,1])-1)*100
  results_vegaf[i, "lci"] <- (exp(confint(tempmod3)[2,1])-1)*100
  results_vegaf[i, "uci"] <- (exp(confint(tempmod3)[2,2])-1)*100
  results_vegaf[i, "p"] <- tempsum3$coefficients[2,4]
  
  # frtaf
  tempmod4 <- lm(paste0(bio[i],"~frtaf_avg+vegaf_avg+ceraf_avg+age_blddraw+race+fast+caco+n2+pmh+asp+mv+nses_avg+
  hbpbase+cholbase+bmi_avg+pa_avg+cigg+pkyrg+alco_avg+calor_avg+poly_avg+sat_avg+trans_avg+gl_avg"), data=tempdata)
  tempsum4 <- summary(tempmod4)
  results_frtaf[i, "n"] <- dim(tempdata)[1]
  results_frtaf[i, "beta"] <- (exp(tempsum4$coefficients[2,1])-1)*100
  results_frtaf[i, "lci"] <- (exp(confint(tempmod4)[2,1])-1)*100
  results_frtaf[i, "uci"] <- (exp(confint(tempmod4)[2,2])-1)*100
  results_frtaf[i, "p"] <- tempsum4$coefficients[2,4]
  
  # stvegaf
  tempmod5 <- lm(paste0(bio[i],"~stvegaf_avg+nonstvegaf_avg+ceraf_avg+frtaf_avg+age_blddraw+race+fast+caco+n2+pmh+asp+mv+nses_avg+
  hbpbase+cholbase+bmi_avg+pa_avg+cigg+pkyrg+alco_avg+calor_avg+poly_avg+sat_avg+trans_avg+gl_avg"), data=tempdata)
  tempsum5 <- summary(tempmod5)
  results_stvegaf[i, "n"] <- dim(tempdata)[1]
  results_stvegaf[i, "beta"] <- (exp(tempsum5$coefficients[2,1])-1)*100
  results_stvegaf[i, "lci"] <- (exp(confint(tempmod5)[2,1])-1)*100
  results_stvegaf[i, "uci"] <- (exp(confint(tempmod5)[2,2])-1)*100
  results_stvegaf[i, "p"] <- tempsum5$coefficients[2,4]

  # nonstveg
  tempmod5 <- lm(paste0(bio[i],"~nonstvegaf_avg+stvegaf_avg+ceraf_avg+frtaf_avg+age_blddraw+race+fast+caco+n2+pmh+asp+mv+nses_avg+
  hbpbase+cholbase+bmi_avg+pa_avg+cigg+pkyrg+alco_avg+calor_avg+poly_avg+sat_avg+trans_avg+gl_avg"), data=tempdata)
  tempsum5 <- summary(tempmod5)
  results_nonstvegaf[i, "n"] <- dim(tempdata)[1]
  results_nonstvegaf[i, "beta"] <- (exp(tempsum5$coefficients[2,1])-1)*100
  results_nonstvegaf[i, "lci"] <- (exp(confint(tempmod5)[2,1])-1)*100
  results_nonstvegaf[i, "uci"] <- (exp(confint(tempmod5)[2,2])-1)*100
  results_nonstvegaf[i, "p"] <- tempsum5$coefficients[2,4]
  
}

results <- rbind(results_aofib,
                 results_ceraf)
results <- rbind(results,
                 results_vegaf)
results <- rbind(results,
                 results_frtaf)
results <- rbind(results,
                 results_stvegaf)
results <- rbind(results,
                 results_nonstvegaf)

results$bon <- results$p*10
results$per <- paste0(format(round(results$beta,digits=1),digits=1),
                      " (",
                      format(round(results$lci,digits=1),digits=1),
                      ", ",
                      format(round(results$uci,digits=1),digits=1),
                      ")")
sum(results$bon<0.05) # 21
View(results[results$bon<0.05,])

save.image(file="/udd/nhyiw/Fiber/biomarker/biomarker.RData")

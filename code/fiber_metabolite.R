library(tidyverse)
library(openxlsx)
metname <- read.xlsx("~/Fiber/metabolomics/metname.xlsx")

load("/udd/nhyiw/Fiber/metabolomics/datapool.RData")
pool=pool2
met275_name <- metname %>% filter(met %in% met275)
f$met <- rownames(f)
met275_name <- merge(f[, c("met", "representative")], met275_name, by="met")
gdata::keep(pool, f, met275_name, met275, sure=T)

#################################################################
#         metabolite-T2D, fiber-metabolite associations
#################################################################
metabolite_t2d <- data.frame(met=met275,
                             name=met275_name$name3,
                             beta=NA,
                             se=NA,
                             HR=NA,
                             p=NA)

aofib_metabolite <- data.frame(met=met275,
                               name=met275_name$name3,
                               beta=NA,
                               se=NA,
                               p=NA)
ceraf_metabolite <- data.frame(met=met275,
                               name=met275_name$name3,
                               beta=NA,
                               se=NA,
                               p=NA)
frtaf_metabolite <- data.frame(met=met275,
                               name=met275_name$name3,
                               beta=NA,
                               se=NA,
                               p=NA)
vegaf_metabolite <- data.frame(met=met275,
                               name=met275_name$name3,
                               beta=NA,
                               se=NA,
                               p=NA)
stvegaf_metabolite <- data.frame(met=met275,
                                 name=met275_name$name3,
                                 beta=NA,
                                 se=NA,
                                 p=NA)
nonstvegaf_metabolite <- data.frame(met=met275,
                                    name=met275_name$name3,
                                    beta=NA,
                                    se=NA,
                                    p=NA)

pool$fast[pool$fast==""] <- "fasting >= 8 hrs"
pool$caco[pool$caco=="not-assig"] <- "control"

pool <- pool %>% mutate(bmi_avg=ifelse(is.na(bmi_avg), median(bmi_avg, na.rm=T), bmi_avg),
                        pa_avg=ifelse(is.na(pa_avg), median(pa_avg, na.rm=T), pa_avg),
                        nses_avg=ifelse(is.na(nses_avg), median(nses_avg, na.rm=T), nses_avg),
                        trans_avg=ifelse(is.na(trans_avg), median(trans_avg, na.rm=T), trans_avg),
                        asp=ifelse(asp==9|is.na(asp), 0, asp),
                        pmh=ifelse(is.na(pmh), 1, pmh),
                        cigg=case_when(is.na(cigg) & pkyrg==1 ~ 1,
                                       is.na(cigg) & pkyrg!=1 ~ 2,
                                       TRUE ~ cigg),
                        aofib_avg=scale(aofib_avg),
                        ceraf_avg=scale(ceraf_avg),
                        frtaf_avg=scale(frtaf_avg),
                        vegaf_avg=scale(vegaf_avg),
                        stvegaf_avg=scale(stvegaf_avg),
                        nonstvegaf_avg=scale(nonstvegaf_avg))

check_na <- data.frame(var=colnames(pool),
                       na=colSums(is.na(pool)))

library(survival)
met_t2d <- list()
met_aofib <- list()
met_fiber <- list()
met_vegfiber <- list()
for (i in 1:length(met275)) {
  tempmet <- met275[i]
  pool[, tempmet] <- scale(pool[, tempmet])
  tempcox <- coxph(as.formula(paste0("Surv(tdb2, type2dbv) ~  ", tempmet, " + age_blddraw + race + fast + nses_avg +
                                     mv + hbpbase + cholbase + pmh + cigg + pkyrg + fmdiab +
                                     alco_avg + bmi_avg + pa_avg + calor_avg + ahei_noal_avg + 
                                     strata(caco, cohort, endpoint)")),
                   data = pool)
  met_t2d[[i]] <- tempcox
  tempsum <- summary(tempcox)
  hr_data <- tempsum$conf.int[1,c("exp(coef)", "lower .95", "upper .95")]
  metabolite_t2d[i, ]$beta <- tempsum$coefficients[1, 1]
  metabolite_t2d[i, ]$se <- tempsum$coefficients[1, 3]
  metabolite_t2d[i, ]$HR <- paste(format(round(hr_data[1], 2), nsmall=2), " (",
                                  format(round(hr_data[2], 2), nsmall=2), ", ",
                                  format(round(hr_data[3], 2), nsmall=2), ")", sep="")
  metabolite_t2d[i, ]$p <- tempsum$coefficients[1, 5]
  
  tempmod <- lm(paste0(tempmet, 
                       "~ aofib_avg + age_blddraw + race + fast + mv + hbpbase + cholbase + 
                       pmh + cigg + pkyrg + nses_avg +
                       alco_avg + bmi_avg + pa_avg + calor_avg + poly_avg + sat_avg + trans_avg +
                       gl_avg + caco + cohort + endpoint"), 
                data=pool)
  met_aofib[[i]] <- tempmod
  aofib_metabolite[i, "beta"] <- summary(tempmod)$coefficients[2,1]
  aofib_metabolite[i, "se"] <- summary(tempmod)$coefficients[2,2]
  aofib_metabolite[i, "p"] <- summary(tempmod)$coefficients[2,4]
  
  tempmod <- lm(paste0(tempmet, 
                       "~ ceraf_avg + frtaf_avg + vegaf_avg + age_blddraw + race + fast + 
                       mv + hbpbase + cholbase + pmh + cigg + pkyrg + nses_avg +
                       alco_avg + bmi_avg + pa_avg + calor_avg + poly_avg + sat_avg + trans_avg +
                       gl_avg + caco + cohort + endpoint"), 
                data=pool)
  met_fiber[[i]] <- tempmod
  ceraf_metabolite[i, "beta"] <- summary(tempmod)$coefficients[2,1]
  ceraf_metabolite[i, "se"] <- summary(tempmod)$coefficients[2,2]
  ceraf_metabolite[i, "p"] <- summary(tempmod)$coefficients[2,4]
  
  frtaf_metabolite[i, "beta"] <- summary(tempmod)$coefficients[3,1]
  frtaf_metabolite[i, "se"] <- summary(tempmod)$coefficients[3,2]
  frtaf_metabolite[i, "p"] <- summary(tempmod)$coefficients[3,4]
  
  vegaf_metabolite[i, "beta"] <- summary(tempmod)$coefficients[4,1]
  vegaf_metabolite[i, "se"] <- summary(tempmod)$coefficients[4,2]
  vegaf_metabolite[i, "p"] <- summary(tempmod)$coefficients[4,4]
  
  tempmod <- lm(paste0(tempmet, 
                       "~ ceraf_avg + frtaf_avg + stvegaf_avg + nonstvegaf_avg + age_blddraw + race + fast + 
                       mv + hbpbase + cholbase + pmh + cigg + pkyrg + nses_avg +
                       alco_avg + bmi_avg + pa_avg + calor_avg + poly_avg + sat_avg + trans_avg +
                       gl_avg + caco + cohort + endpoint"), 
                data=pool)
  met_vegfiber[[i]] <- tempmod
  stvegaf_metabolite[i, "beta"] <- summary(tempmod)$coefficients[4,1]
  stvegaf_metabolite[i, "se"] <- summary(tempmod)$coefficients[4,2]
  stvegaf_metabolite[i, "p"] <- summary(tempmod)$coefficients[4,4]
  
  nonstvegaf_metabolite[i, "beta"] <- summary(tempmod)$coefficients[5,1]
  nonstvegaf_metabolite[i, "se"] <- summary(tempmod)$coefficients[5,2]
  nonstvegaf_metabolite[i, "p"] <- summary(tempmod)$coefficients[5,4]
}

########################################################
#                   outcomes
########################################################
metabolite_t2d$fdr <- p.adjust(metabolite_t2d$p, method="fdr", length(metabolite_t2d$p))
metabolite_t2d$bon <- p.adjust(metabolite_t2d$p, method="bonferroni", length(metabolite_t2d$p))
sum(metabolite_t2d$bon<0.05) # 120

length(Reduce(union, list(aofib_metabolite %>% mutate(bon=p.adjust(p, method="bonferroni", n=275)) %>% 
                            filter(bon<0.05) %>% pull(met),
                          ceraf_metabolite %>% mutate(bon=p.adjust(p, method="bonferroni", n=275)) %>% 
                            filter(bon<0.05) %>% pull(met),
                          frtaf_metabolite %>% mutate(bon=p.adjust(p, method="bonferroni", n=275)) %>% 
                            filter(bon<0.05) %>% pull(met),
                          vegaf_metabolite %>% mutate(bon=p.adjust(p, method="bonferroni", n=275)) %>% 
                            filter(bon<0.05) %>% pull(met),
                          nonstvegaf_metabolite %>% mutate(bon=p.adjust(p, method="bonferroni", n=275)) %>% 
                            filter(bon<0.05) %>% pull(met),
                          stvegaf_metabolite %>% mutate(bon=p.adjust(p, method="bonferroni", n=275)) %>% 
                            filter(bon<0.05) %>% pull(met)))) # 120

aofib_metabolite <- aofib_metabolite %>% 
  mutate(bon=p.adjust(p, method="bonferroni", n=275))
ceraf_metabolite <- ceraf_metabolite %>% 
  mutate(bon=p.adjust(p, method="bonferroni", n=275))
frtaf_metabolite <- frtaf_metabolite %>% 
  mutate(bon=p.adjust(p, method="bonferroni", n=275))
vegaf_metabolite <- vegaf_metabolite %>% 
  mutate(bon=p.adjust(p, method="bonferroni", n=275))
stvegaf_metabolite <- stvegaf_metabolite %>% 
  mutate(bon=p.adjust(p, method="bonferroni", n=275))
nonstvegaf_metabolite <- nonstvegaf_metabolite %>% 
  mutate(bon=p.adjust(p, method="bonferroni", n=275))

gdata::keep(met275_name, metabolite_t2d, aofib_metabolite, ceraf_metabolite, frtaf_metabolite,
            vegaf_metabolite, stvegaf_metabolite, nonstvegaf_metabolite, sure=T)
save.image(file="/udd/nhyiw/Fiber/metabolomics/t2d/fiber_metabolite.RData")
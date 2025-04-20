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
pool <- pool %>% filter(canbase==0 & hrtbase==0 & strbase==0 & dbbase==0 & is.na(dbbase2) & caco==0)

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
  tempmod1 <- lm(paste0(bio[i],"~aofib_avg+age_blddraw+race+fast+n2+pmh+asp+mv+hbpbase+cholbase+nses_avg+
  bmi_avg+pa_avg+cigg+pkyrg+alco_avg+calor_avg+poly_avg+sat_avg+trans_avg+gl_avg"), data=tempdata)
  tempsum1 <- summary(tempmod1)
  results_aofib[i, "n"] <- dim(tempdata)[1]
  results_aofib[i, "beta"] <- (exp(tempsum1$coefficients[2,1])-1)*100
  results_aofib[i, "lci"] <- (exp(confint(tempmod1)[2,1])-1)*100
  results_aofib[i, "uci"] <- (exp(confint(tempmod1)[2,2])-1)*100
  results_aofib[i, "p"] <- tempsum1$coefficients[2,4]
  
  # ceraf
  tempmod2 <- lm(paste0(bio[i],"~ceraf_avg+frtaf_avg+vegaf_avg+age_blddraw+race+fast+n2+pmh+asp+mv+nses_avg+
  hbpbase+cholbase+bmi_avg+pa_avg+cigg+pkyrg+alco_avg+calor_avg+poly_avg+sat_avg+trans_avg+gl_avg"), data=tempdata)
  tempsum2 <- summary(tempmod2)
  results_ceraf[i, "n"] <- dim(tempdata)[1]
  results_ceraf[i, "beta"] <- (exp(tempsum2$coefficients[2,1])-1)*100
  results_ceraf[i, "lci"] <- (exp(confint(tempmod2)[2,1])-1)*100
  results_ceraf[i, "uci"] <- (exp(confint(tempmod2)[2,2])-1)*100
  results_ceraf[i, "p"] <- tempsum2$coefficients[2,4]
  
  # vegaf
  tempmod3 <- lm(paste0(bio[i],"~vegaf_avg+ceraf_avg+frtaf_avg+age_blddraw+race+fast+n2+pmh+asp+mv+nses_avg+
  hbpbase+cholbase+bmi_avg+pa_avg+cigg+pkyrg+alco_avg+calor_avg+poly_avg+sat_avg+trans_avg+gl_avg"), data=tempdata)
  tempsum3 <- summary(tempmod3)
  results_vegaf[i, "n"] <- dim(tempdata)[1]
  results_vegaf[i, "beta"] <- (exp(tempsum3$coefficients[2,1])-1)*100
  results_vegaf[i, "lci"] <- (exp(confint(tempmod3)[2,1])-1)*100
  results_vegaf[i, "uci"] <- (exp(confint(tempmod3)[2,2])-1)*100
  results_vegaf[i, "p"] <- tempsum3$coefficients[2,4]
  
  # frtaf
  tempmod4 <- lm(paste0(bio[i],"~frtaf_avg+vegaf_avg+ceraf_avg+age_blddraw+race+fast+n2+pmh+asp+mv+nses_avg+
  hbpbase+cholbase+bmi_avg+pa_avg+cigg+pkyrg+alco_avg+calor_avg+poly_avg+sat_avg+trans_avg+gl_avg"), data=tempdata)
  tempsum4 <- summary(tempmod4)
  results_frtaf[i, "n"] <- dim(tempdata)[1]
  results_frtaf[i, "beta"] <- (exp(tempsum4$coefficients[2,1])-1)*100
  results_frtaf[i, "lci"] <- (exp(confint(tempmod4)[2,1])-1)*100
  results_frtaf[i, "uci"] <- (exp(confint(tempmod4)[2,2])-1)*100
  results_frtaf[i, "p"] <- tempsum4$coefficients[2,4]
  
  # stvegaf
  tempmod5 <- lm(paste0(bio[i],"~stvegaf_avg+nonstvegaf_avg+ceraf_avg+frtaf_avg+age_blddraw+race+fast+n2+pmh+asp+mv+nses_avg+
  hbpbase+cholbase+bmi_avg+pa_avg+cigg+pkyrg+alco_avg+calor_avg+poly_avg+sat_avg+trans_avg+gl_avg"), data=tempdata)
  tempsum5 <- summary(tempmod5)
  results_stvegaf[i, "n"] <- dim(tempdata)[1]
  results_stvegaf[i, "beta"] <- (exp(tempsum5$coefficients[2,1])-1)*100
  results_stvegaf[i, "lci"] <- (exp(confint(tempmod5)[2,1])-1)*100
  results_stvegaf[i, "uci"] <- (exp(confint(tempmod5)[2,2])-1)*100
  results_stvegaf[i, "p"] <- tempsum5$coefficients[2,4]
  
  # nonstveg
  tempmod5 <- lm(paste0(bio[i],"~nonstvegaf_avg+stvegaf_avg+ceraf_avg+frtaf_avg+age_blddraw+race+fast+n2+pmh+asp+mv+nses_avg+
  hbpbase+cholbase+bmi_avg+pa_avg+cigg+pkyrg+alco_avg+calor_avg+poly_avg+sat_avg+trans_avg+gl_avg"), data=tempdata)
  tempsum5 <- summary(tempmod5)
  results_nonstvegaf[i, "n"] <- dim(tempdata)[1]
  results_nonstvegaf[i, "beta"] <- (exp(tempsum5$coefficients[2,1])-1)*100
  results_nonstvegaf[i, "lci"] <- (exp(confint(tempmod5)[2,1])-1)*100
  results_nonstvegaf[i, "uci"] <- (exp(confint(tempmod5)[2,2])-1)*100
  results_nonstvegaf[i, "p"] <- tempsum5$coefficients[2,4]
  
}

results_con <- rbind(results_aofib,
                 results_ceraf)
results_con <- rbind(results_con,
                 results_vegaf)
results_con <- rbind(results_con,
                 results_frtaf)
results_con <- rbind(results_con,
                 results_stvegaf)
results_con <- rbind(results_con,
                 results_nonstvegaf)

results_con$bon <- results_con$p*10
results_con$per <- paste0(format(round(results_con$beta,digits=1),digits=1),
                          " (",
                          format(round(results_con$lci,digits=1),digits=1),
                          ", ",
                          format(round(results_con$uci,digits=1),digits=1),
                          ")")
sum(results_con$bon<0.05) # 14
View(results_con[results_con$p<0.05,])
View(results_con[results_con$p<0.005,])

save.image(file="/udd/nhyiw/Fiber/biomarker/biomarker_con.RData")

load("~/Fiber/biomarker/biomarker_con.RData")
gdata::keep(results_con, sure=T)
load("~/Fiber/biomarker/biomarker.RData")
gdata::keep(results_con, results, sure=T)

results <- results %>% mutate(bio=case_when(bio=="adj_lna1c" ~ "HbA1c",
                                            bio=="adj_lnadipo" ~ "Adiponectin",
                                            bio=="adj_lnchol" ~ "Total cholesterol",
                                            bio=="adj_lncpep" ~ "C-peptide",
                                            bio=="adj_lncrp" ~ "CRP",
                                            bio=="adj_lnhdl" ~ "HDL-C",
                                            bio=="adj_lnldl" ~ "LDL-C",
                                            bio=="adj_lnleptin" ~ "Leptin",
                                            bio=="adj_lnratio" ~ "TAG/HDL-C ratio",
                                            bio=="adj_lntrig" ~ "TAG"))
results_con <- results_con %>% mutate(bio=case_when(bio=="adj_lna1c" ~ "HbA1c",
                                            bio=="adj_lnadipo" ~ "Adiponectin",
                                            bio=="adj_lnchol" ~ "Total cholesterol",
                                            bio=="adj_lncpep" ~ "C-peptide",
                                            bio=="adj_lncrp" ~ "CRP",
                                            bio=="adj_lnhdl" ~ "HDL-C",
                                            bio=="adj_lnldl" ~ "LDL-C",
                                            bio=="adj_lnleptin" ~ "Leptin",
                                            bio=="adj_lnratio" ~ "TAG/HDL-C ratio",
                                            bio=="adj_lntrig" ~ "TAG"))

aofib0 <- results %>% filter(type=="aofib")
aofib1 <- results_con %>% filter(type=="aofib")
ceraf0 <- results %>% filter(type=="ceraf")
ceraf1 <- results_con %>% filter(type=="ceraf")
frtaf0 <- results %>% filter(type=="frtaf")
frtaf1 <- results_con %>% filter(type=="frtaf")
stvegaf0 <- results %>% filter(type=="stvegaf")
stvegaf1 <- results_con %>% filter(type=="stvegaf")

aofib <- data.frame(type="aofib", bio=aofib0$bio, beta1=aofib0$beta, beta2=aofib1$beta, 
                    group=ifelse(aofib0$bon<0.05, 
                                 "Significant (Bonferroni-adjusted P <0.05)", 
                                 "Nonsignificant")) %>% 
  mutate(group=factor(group, levels=c("Significant (Bonferroni-adjusted P <0.05)", 
                                      "Nonsignificant")))
ceraf <- data.frame(type="ceraf", bio=ceraf0$bio, beta1=ceraf0$beta, beta2=ceraf1$beta, 
                    group=ifelse(ceraf0$bon<0.05, 
                                 "Significant (Bonferroni-adjusted P <0.05)", 
                                 "Nonsignificant")) %>% 
  mutate(group=factor(group, levels=c("Significant (Bonferroni-adjusted P <0.05)", 
                                      "Nonsignificant")))
frtaf <- data.frame(type="frtaf", bio=frtaf0$bio, beta1=frtaf0$beta, beta2=frtaf1$beta, 
                    group=ifelse(frtaf0$bon<0.05, 
                                 "Significant (Bonferroni-adjusted P <0.05)", 
                                 "Nonsignificant")) %>% 
  mutate(group=factor(group, levels=c("Significant (Bonferroni-adjusted P <0.05)", 
                                      "Nonsignificant")))
stvegaf <- data.frame(type="stvegaf", bio=stvegaf0$bio, beta1=stvegaf0$beta, beta2=stvegaf1$beta, 
                      group=ifelse(stvegaf0$bon<0.05, 
                                   "Significant (Bonferroni-adjusted P <0.05)", 
                                   "Nonsignificant")) %>% 
  mutate(group=factor(group, levels=c("Significant (Bonferroni-adjusted P <0.05)", 
                                      "Nonsignificant")))

figS1 <- bind_rows(list(aofib, ceraf, frtaf, stvegaf))
dotcolor <- c("#d7191c", "grey")
names(dotcolor) <- c("Significant (Bonferroni-adjusted P <0.05)", 
                     "Nonsignificant")

library(ggrepel)
fb <- c("aofib", "ceraf", "frtaf", "stvegaf")
fiber <- c("Total fiber", "Fiber from cereals", "Fiber from fruits",
           "Fiber from starchy vegetables")
pp <- list()
for (i in 1:4) {
  temp=figS1 %>% filter(type==fb[i])
  
  pp[[i]]=ggplot(aes(x=beta1, y=beta2, label=bio), data=temp) +
    geom_abline(slope=1, linetype=2, size=1.5) +
    geom_hline(yintercept = 0, linetype=2, size=1.5) +
    geom_vline(xintercept = 0, linetype=2, size=1.5) +
    geom_point(aes(color=group), size=3.5) +
    scale_y_continuous(breaks=c(-10, -5, 0, 5, 10), limits=c(-12, 10)) +
    scale_x_continuous(breaks=c(-10, -5, 0, 5, 10), limits=c(-12, 10)) +
    scale_color_manual(values=dotcolor) +
    labs(
      title=fiber[i],
      x="% Difference\namong all participants",
      y="% Difference\namong controls") +
    theme_bw() +
    theme(panel.grid = element_blank(),
          panel.border = element_rect(size=2.2)) +
    theme(plot.title = element_text(hjust = 0, size=18, color="black"),
          axis.title = element_text(size=18, color="black"),
          axis.text = element_text(size=18, color="black"),
          axis.ticks = element_line(size=1.2),
          axis.ticks.length = unit(.2, "cm"),
          legend.position = "none") +
    geom_label_repel(size=5,
                     box.padding   = 0.35, 
                     point.padding = 0.5,
                     segment.color = 'grey50',
                     max.overlaps = Inf) +
    theme(aspect.ratio = 1)
}

pdf("/udd/nhyiw/Fiber/biomarker/figS1.pdf", 
    width = 12, height = 12, onefile = F) # Open a new pdf file
egg::ggarrange(pp[[1]], pp[[2]], pp[[3]], pp[[4]],
               nrow = 2)
dev.off() # Close the file

figS1_legend=ggplot(aes(x=beta1, y=beta2, label=bio), data=temp) +
  geom_abline(slope=1, linetype=2, size=1.5) +
  geom_hline(yintercept = 0, linetype=2, size=1.5) +
  geom_vline(xintercept = 0, linetype=2, size=1.5) +
  geom_point(aes(color=group), size=3.5) +
  scale_y_continuous(breaks=c(-10, -5, 0, 5, 10), limits=c(-12, 10)) +
  scale_x_continuous(breaks=c(-10, -5, 0, 5, 10), limits=c(-12, 10)) +
  scale_color_manual(values=dotcolor) +
  labs(
    title=fiber[i],
    x="% Difference\namong all participants",
    y="% Difference\namong controls",
    color="Significance of heme iron-biomarker \nassociations among all participants") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        panel.border = element_rect(size=2.2)) +
  theme(plot.title = element_text(hjust = 0, size=18, color="black"),
        axis.title = element_text(size=18, color="black"),
        axis.text = element_text(size=18, color="black"),
        axis.ticks = element_line(size=1.2),
        axis.ticks.length = unit(.2, "cm"),
        legend.title = element_text(size = 16),
        legend.text = element_text(size = 16),
        legend.position = "right")

pdf("/udd/nhyiw/Fiber/biomarker/figS1_legend.pdf", 
    width = 6, height = 3, onefile = F) # Open a new pdf file
grid::grid.draw(ggpubr::get_legend(figS1_legend))
dev.off() # Close the file

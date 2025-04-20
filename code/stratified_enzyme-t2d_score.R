# only among those with both mgx and mbx
library(tidyverse)
library(lmerTest)
setwd("/n/home03/ywan/fiber/0430")

load("/n/home03/ywan/fiber/_data/pool_enzyme_spe.RData")

##############################################################################
#                    run lm for the t2dscore
##############################################################################
ecs_data <- data.frame(t(pool_enzyme))
identical(rownames(ecs_data), metadata_INT$id) # TRUE

dataall0 <- cbind(metadata_INT[, c(1:7,11:16,19:34)], ecs_data)
dataall2 <- dataall0[, -2]
dataall2 <- dataall2 %>% group_by(randid, study, sex) %>% summarise_all(mean, na.rm=T)

load("~/fiber/0430/t2dscore.RData")
dataall2 <- merge(dataall2, metadata_all[, c("randid", "t2dscore", "t2dscore_std")], 
                  by="randid")

# log transform species and biomarker
LOG <- function(x) {
  y <- replace(x, x == 0, min(x[x>0]) / 2)
  return(log10(y))
}
for (i in 29:4648){
  dataall2[, i] <- LOG(dataall2[, i]*100)
}

data_mbs2 <- dataall2[dataall2$study=="NHSII (USA)", ]
data_mlvs2 <- dataall2[dataall2$study=="HPFS (USA)", ]
for (i in 29:4648){
  data_mbs2[, i] <- scale(data_mbs2[, i])
}
for (i in 29:4648){
  data_mlvs2[, i] <- scale(data_mlvs2[, i])
}

score_res_mbs <- data.frame(bio="score",
                            pair=colnames(dataall2)[29:4648],
                            beta=NA,
                            se=NA,
                            p=NA)
score_res_mlvs <- score_res_mbs

for (j in 1:4620){
  mbs_mod <- summary(lm(as.formula(paste0("t2dscore_std~", 
                                          colnames(dataall2)[j+28], 
                                          "+age+bmi+antib+ahei_noal")),
                        data=data_mbs2))
  mlvs_mod <- summary(lm(as.formula(paste0("t2dscore_std~", 
                                           colnames(dataall2)[j+28], 
                                           "+age+bmi+antib+ahei_noal")),
                         data=data_mlvs2))
  
  score_res_mbs[j, "beta"] <- mbs_mod$coefficients[2,1]
  score_res_mbs[j, "se"] <- mbs_mod$coefficients[2,2]  
  score_res_mbs[j, "p"] <- mbs_mod$coefficients[2,4]
  
  score_res_mlvs[j, "beta"] <- mlvs_mod$coefficients[2,1]
  score_res_mlvs[j, "se"] <- mlvs_mod$coefficients[2,2]  
  score_res_mlvs[j, "p"] <- mlvs_mod$coefficients[2,4]
}

score_res_mbs$fdr <- p.adjust(score_res_mbs$p, method="fdr", n=length(score_res_mbs$p))
score_res_mlvs$fdr <- p.adjust(score_res_mlvs$p, method="fdr", n=length(score_res_mlvs$p))

score_res_mlvs <- merge(enzyme, score_res_mlvs, by="pair")
score_res_mbs <- merge(enzyme, score_res_mbs, by="pair")

##############################################################################
#                    run meta-analysis
##############################################################################
gdata::keep(score_res_mlvs, score_res_mbs, sure=T)
score_res_mbs$study <- "MBS"
score_res_mlvs$study <- "MLVS"

meta_bio_res <- score_res_mlvs[,1:5]

library(meta)
for(i in 1:dim(score_res_mlvs)[1]){
  temp <- rbind(score_res_mlvs[i,],
                score_res_mbs[i,])
  
  dd <- metagen(beta,
                se,
                studlab = study,
                method.tau = "DL",
                sm = "MD",
                data = temp)
  meta_bio_res[i, "beta"] <- dd$TE.fixed
  meta_bio_res[i, "se"] <- dd$seTE.fixed
  meta_bio_res[i, "p"] <- dd$pval.fixed
}

meta_bio_res$fdr <- p.adjust(meta_bio_res$p, method="fdr", n=4620)
View(meta_bio_res[meta_bio_res$fdr<0.25,])
table(meta_bio_res[meta_bio_res$fdr<0.25,]$species)

save.image(file="/n/home03/ywan/fiber/0430/ecs_bio_mv.RData")

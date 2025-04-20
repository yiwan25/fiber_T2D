# only among those with both mgx and mbx
library(tidyverse)
library(lmerTest)
setwd("/n/home03/ywan/fiber/0430")

load(file="/n/home03/ywan/fiber/0430/fiber_spe.RData")

##############################################################################
#                 run lmer for the biomarker
##############################################################################
spe_data <- data.frame(t(species_zero))
identical(rownames(spe_data), metadata_INT$id) # TRUE

dataall0 <- cbind(metadata_INT[, c(1:7,19:34)], spe_data)
dataall <- dataall0

# log transform species and biomarker
LOG <- function(x) {
  y <- replace(x, x == 0, min(x[x>0]) / 2)
  return(log10(y))
}

##############################################################################
#                    run lm for the t2dscore
##############################################################################
dataall2 <- dataall0[, -2]
dataall2 <- dataall2 %>% group_by(randid, study, sex) %>% summarise_all(mean, na.rm=T)

load("~/fiber/0430/t2dscore.RData")
dataall2 <- merge(dataall2, metadata_all[, c("randid", "t2dscore", "t2dscore_std")], 
                  by="randid")

# log transform species and biomarker
for (i in 1:139){
  dataall2[, spe_final[i]] <- LOG(dataall2[, spe_final[i]])
}

data_mbs2 <- dataall2[dataall2$study=="NHSII (USA)", ]
data_mlvs2 <- dataall2[dataall2$study=="HPFS (USA)", ]
for (i in 1:139){
  data_mbs2[, spe_final[i]] <- scale(data_mbs2[, spe_final[i]])
}
for (i in 1:139){
  data_mlvs2[, spe_final[i]] <- scale(data_mlvs2[, spe_final[i]])
}

score_res_mbs <- data.frame(bio="score",
                            spe=spe_final,
                            beta=NA,
                            se=NA,
                            p=NA)
score_res_mlvs <- score_res_mbs

for (j in 1:length(spe_final)){
  mbs_mod <- summary(lm(as.formula(paste0("t2dscore_std~", 
                                          spe_final[j], 
                                          "+age+bmi+antib+ahei_noal")),
                        data=data_mbs2))
  mlvs_mod <- summary(lm(as.formula(paste0("t2dscore_std~", 
                                           spe_final[j], 
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

save.image(file="/n/home03/ywan/fiber/0430/spe_score_mv.RData")

# figure 4b
spe <- c("s__Clostridium_bolteae", "s__Clostridium_symbiosum", "s__Ruminococcus_gnavus")
spelab <- c("Clostridium bolteae", "Clostridium symbiosum", "Ruminococcus gnavus")

library(ggplot2)
pp2 <- list()
for (i in 1:3) {
  aa <- data.frame(x=dataall2[, spe[i]],
                   y=dataall2$t2dscore_std)
  # aa2 <- aa[aa$spe>min(aa$spe), ]
  
  pp2[[i]] <-   ggplot(aes(x=x, y=y), data=aa) +
    geom_point(fill="#d01c8b", shape=21, alpha=0.3) +
    geom_smooth(col="black", method="lm", 
                # aes(x=fib, y=spe), data=aa2,
                se=F) +
    labs(x="Log10 relative abundance",
         y="Standardized T2D metabolomic score",
         title = "") +
    theme_bw() +
    theme(panel.grid = element_blank(),
          panel.border = element_rect(linewidth = 1),
          axis.title = element_text(colour = "black", size=15),
          axis.text = element_text(colour = "black", size=15),
          plot.title = element_text(colour = "black", size=15, hjust = 0.5,
                                    face = "italic"))
}

load(file="/n/home03/ywan/fiber/0430/fiber_spe_fig4a.RData")

pdf("/n/home03/ywan/fiber/0430/fig4a.pdf", 
    width = 11, height = 8, onefile = F) # Open a new pdf file
egg::ggarrange(pp[[1]], pp[[2]], pp[[3]], 
               pp2[[1]], pp2[[2]], pp2[[3]], nrow = 2)
dev.off() # Close the file

##############################################################################
#                    run meta-analysis
##############################################################################
gdata::keep(score_res_mbs, score_res_mlvs, sure=T)
score_res_mbs$study <- "MBS"
score_res_mlvs$study <- "MLVS"

identical(paste0(score_res_mbs$bio, "_", score_res_mbs$spe),
          paste0(score_res_mlvs$bio, "_", score_res_mlvs$spe)) # TRUE

meta_score_res <- score_res_mlvs[,1:5]

library(meta)
for(i in 1:dim(score_res_mbs)[1]){
  temp <- rbind(score_res_mbs[i,],
                score_res_mlvs[i,])
  
  dd <- metagen(beta,
                se,
                studlab = study,
                method.tau = "DL",
                sm = "MD",
                data = temp)
  meta_score_res[i, "beta"] <- dd$TE.fixed
  meta_score_res[i, "se"] <- dd$seTE.fixed
  meta_score_res[i, "p"] <- dd$pval.fixed
}

meta_score_res2 <- data.frame()
name <- unique(meta_score_res$bio)
for (i in 1:length(name)){
  temp <- meta_score_res[meta_score_res$bio==name[i],]
  temp$fdr <- p.adjust(temp$p, method="fdr", n=length(temp$p))
  
  meta_score_res2 <- rbind(meta_score_res2, temp)
}
save.image(file="/n/home03/ywan/fiber/0430/spe_score_mv2.RData")

# only among those with both mgx and mbx
library(tidyverse)
library(MMUPHin)
setwd("/n/home03/ywan/fiber/0430")

load("/n/home03/ywan/fiber/_data/pool_species_dt.RData")

# meta-analysis
species_zero <- species_zero/100

# can't use _zero to calculate relab
spe_relab <- data.frame(species=rownames(species),
                        relab=rowMeans(species, na.rm = T)) %>% arrange(desc(relab))

metadata$study <- as.factor(metadata$study)
identical(metadata$id, colnames(species_zero))
rownames(metadata) <- metadata$id

table(spe_all$present)
# 1   2 
# 24 139

# exclude solo species
spe_final <- spe_all[spe_all$present>1,]$species
species_zero0 <- species_zero
species_zero <- species_zero[(rownames(species_zero) %in% spe_final),]

summary(colSums(species_zero, na.rm = T))
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.4683  0.9388  0.9787  0.9521  0.9939  1.0000

# add "other" category
species_zero[140,] <- t(1-colSums(species_zero, na.rm = T))
rownames(species_zero)[140] <- "s__Other"
summary(colSums(species_zero, na.rm = T)) # should be all 1
sum(species_zero<0, na.rm = T) # 7
species_zero[species_zero<0] <- 0

summary(metadata$age)
summary(metadata$bmi)
table(metadata$sex)
metadata$study <- as.factor(metadata$study)
metadata$bmi <- as.numeric(metadata$bmi)

fbgp <- c("aofib", "ceraf", "frtaf", "vegaf", "nonstvegaf", "stvegaf")

# INT within each study
temp <- metadata[!duplicated(metadata$randid), ]
temp_INT <- temp

library(gtools)
INT_trans = function(aa) {
  bb <- qnorm((rank(aa, na.last="keep") - 0.5) / sum(!is.na(aa)))
  return(bb)
}
for (i in 1:length(fbgp)) {
  # the ave function makes the INT_trans function stratifed by study
  temp_INT[, fbgp[i]] <- ave(temp_INT[, fbgp[i]], temp_INT$study, FUN = INT_trans)
}

metadata_INT <- merge(metadata[, c(1:14,49,15:18,25:33)], temp_INT[, c("randid", fbgp)], by="randid")
species_zero <- species_zero[, metadata_INT$id]
rownames(metadata_INT) <- metadata_INT$id

save.image(file="/n/home03/ywan/fiber/0430/fiber_spe.RData")

##############################################################################
#                 run lm_meta for the main dataset
##############################################################################
covar <- list()
covar[[1]] <- c("calor")
covar[[2]] <- c("calor", "frtaf", "vegaf")
covar[[3]] <- c("calor", "ceraf", "vegaf")
covar[[4]] <- c("calor", "ceraf", "frtaf")
covar[[5]] <- c("calor", "ceraf", "frtaf", "stvegaf")
covar[[6]] <- c("calor", "ceraf", "frtaf", "nonstvegaf")

species_zero[is.na(species_zero)] <- 0
identical(colnames(species_zero), rownames(metadata_INT)) # TRUE
meta_fiber_INT <- list()
meta_results_INT <- list()

for (i in 1:length(fbgp)) {
  meta_fiber_INT[[i]] <- lm_meta(feature_abd = as.matrix(species_zero),
                                 batch = "study",
                                 exposure = fbgp[i],
                                 covariates = c("age", "bmi", "antib", "ahei2010", covar[[i]]),
                                 covariates_random = 'randid',
                                 data = metadata_INT,
                                 control = list(verbose = T,
                                                rma_method = "FE",
                                                transform = "LOG",
                                                output = paste0("./fiber_spe_mv/", fbgp[i])))
  meta_results_INT[[i]] <- meta_fiber_INT[[i]]$meta_fits
}

identical(colnames(meta_results_INT[[1]]), colnames(meta_results_INT[[2]])) # TRUE
speresults_INT <- bind_rows(meta_results_INT)
speresults_INT$feature <- gsub("s__", "", speresults_INT$feature)

table(speresults_INT[speresults_INT$qval.fdr<0.25,]$exposure)

save.image(file="/n/home03/ywan/fiber/0430/fiber_spe_mv.RData")

# figure 3
spe <- c("s__Bifidobacterium_adolescentis", "s__Eubacterium_eligens", "s__Gemmiger_formicilis")
spelab <- c("Bifidobacterium adolescentis", "Eubacterium eligens", "Gemmiger formicilis")
fib <- c("ceraf", "frtaf", "frtaf")
fibcol <- c("#d7191c", "#fdae61", "#fdae61")
fiblab <- c("Fiber from cereals", 
            "Fiber from fruits", 
            "Fiber from fruits")

identical(colnames(species_zero), metadata$id) # TRUE
identical(colnames(species_zero), metadata_INT$id) # TRUE

LOG <- function(x) {
  y <- replace(x, x == 0, min(x[x>0]) / 2)
  return(log10(y))
}

library(ggplot2)
pp <- list()
for (i in 1:3) {
  aa <- data.frame(x=metadata[, fib[i]],
                   y=t(species_zero[rownames(species_zero)==spe[i], ]))
  colnames(aa) <- c("fib", "spe")
  aa$spe <- LOG(aa$spe)
  # aa2 <- aa[aa$spe>min(aa$spe), ]
  
  pp[[i]] <-   ggplot(aes(x=fib, y=spe), data=aa) +
    geom_point(fill=fibcol[i], shape=21, alpha=0.3) +
    geom_smooth(col="black", method="lm", 
                # aes(x=fib, y=spe), data=aa2,
                se=F) +
    labs(x=fiblab[i],
         y="Log10 relative abundance",
         title = spelab[i]) +
    theme_bw() +
    theme(panel.grid = element_blank(),
          panel.border = element_rect(linewidth = 1),
          axis.title = element_text(colour = "black", size=15),
          axis.text = element_text(colour = "black", size=15),
          plot.title = element_text(colour = "black", size=15, hjust = 0.5,
                                    face = "italic"))
}

pdf("/n/home03/ywan/fiber/0430/fiber_spe.pdf", 
    width = 4, height = 12, onefile = F) # Open a new pdf file
egg::ggarrange(pp[[1]], pp[[2]], pp[[3]], nrow = 3)
dev.off() # Close the file

# figure 3c
spe <- c("s__Clostridium_bolteae", "s__Clostridium_symbiosum", "s__Ruminococcus_gnavus")
spelab <- c("Clostridium bolteae", "Clostridium symbiosum", "Ruminococcus gnavus")
fib <- c("frtaf", "frtaf", "frtaf")
fibcol <- c("#fdae61", "#fdae61", "#fdae61")
fiblab <- c("Fiber from fruits", 
            "Fiber from fruits", 
            "Fiber from fruits")

identical(colnames(species_zero), metadata$id) # TRUE
identical(colnames(species_zero), metadata_INT$id) # TRUE

LOG <- function(x) {
  y <- replace(x, x == 0, min(x[x>0]) / 2)
  return(log10(y))
}

library(ggplot2)
pp <- list()
for (i in 1:3) {
  aa <- data.frame(x=metadata[, fib[i]],
                   y=t(species_zero[rownames(species_zero)==spe[i], ]))
  colnames(aa) <- c("fib", "spe")
  aa$spe <- LOG(aa$spe)
  # aa2 <- aa[aa$spe>min(aa$spe), ]
  
  pp[[i]] <-   ggplot(aes(x=fib, y=spe), data=aa) +
    geom_point(fill=fibcol[i], shape=21, alpha=0.3) +
    geom_smooth(col="black", method="lm", 
                # aes(x=fib, y=spe), data=aa2,
                se=F) +
    labs(x=fiblab[i],
         y="Log10 relative abundance",
         title = spelab[i]) +
    theme_bw() +
    theme(panel.grid = element_blank(),
          panel.border = element_rect(linewidth = 1),
          axis.title = element_text(colour = "black", size=15),
          axis.text = element_text(colour = "black", size=15),
          plot.title = element_text(colour = "black", size=15, hjust = 0.5,
                                    face = "italic"))
}

save(pp, file="/n/home03/ywan/fiber/0430/fiber_spe_fig4a.RData")

pdf("/n/home03/ywan/fiber/0430/fiber_spe_fig3c.pdf", 
    width = 13, height = 4, onefile = F) # Open a new pdf file
egg::ggarrange(pp[[1]], pp[[2]], pp[[3]], nrow = 1)
dev.off() # Close the file

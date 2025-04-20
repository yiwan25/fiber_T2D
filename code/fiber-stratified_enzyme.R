# only among those with both mgx and mbx
library(tidyverse)
library(MMUPHin)

setwd("/n/home03/ywan/fiber/0430")
load("/n/home03/ywan/fiber/_data/pool_enzyme_spe.RData")
fbgp <- c("aofib", "ceraf", "frtaf", "vegaf", "nonstvegaf", "stvegaf")

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

meta_fiber_INT <- list()
meta_results_INT <- list()

library(future)
future::plan(multisession, workers = 86)
for (i in 1:length(fbgp)) {
  meta_fiber_INT[[i]] <- lm_meta(feature_abd = as.matrix(pool_enzyme),
                                 batch = "study",
                                 exposure = fbgp[i],
                                 covariates = c("age", "bmi", "antib", "ahei2010", covar[[i]]),
                                 covariates_random = 'randid',
                                 data = metadata_INT,
                                 control = list(verbose = T,
                                                rma_method = "FE",
                                                transform = "LOG",
                                                normalization = "None",
                                                output = paste0("./fiber_ecs_mv/", fbgp[i])))
  meta_results_INT[[i]] <- meta_fiber_INT[[i]]$meta_fits
}

identical(colnames(meta_results_INT[[1]]), colnames(meta_results_INT[[2]])) # TRUE
ecsresults_INT <- bind_rows(meta_results_INT)
ecsresults_INT$feature <- gsub("s__", "", ecsresults_INT$feature)
rownames(ecsresults_INT) <- NULL
colnames(ecsresults_INT)[1] <- "pair"

enzy <- data.frame(str_split(enzyme$enzy, "\\:", simplify = T))
enzyme$ec <- enzy$X1
ecsresults_INT <- merge(enzyme, ecsresults_INT, by="pair")

save.image(file="/n/home03/ywan/fiber/0430/fiber_ecs_mv.RData")

length(unique(ecsresults_INT$species)) # 16 out of 21 (5 do not have ec)
View(ecsresults_INT[ecsresults_INT$qval.fdr<0.25,])
table(ecsresults_INT[ecsresults_INT$qval.fdr<0.25 & ecsresults_INT$exposure=="frtaf",]$species)
table(ecsresults_INT[ecsresults_INT$ec %in% c(pec_ec, cel_ec) & 
                       ecsresults_INT$qval.fdr<0.25,]$species)

# pectin degradation
pec_ec <- c("EC3.2.1.23", "EC3.2.1.55",
            "EC3.1.1.11", "EC3.2.1.82", "EC3.2.1.89", "EC3.2.1.172", "EC4.2.2.2", "EC5.3.1.12")

# cellulose degradation
cel_ec <- c("EC2.4.1.281", "EC3.2.1.4", "EC3.2.1.8", "EC3.2.1.21", 
            "EC3.2.1.37", "EC3.2.1.40", "EC5.1.3.11",
            "EC3.2.1.23",
            "EC3.2.1.55")

ec_all <- data.frame(ec=c("EC3.2.1.40", "EC3.2.1.4", "EC5.1.3.11", "EC2.4.1.281",
                          "EC3.2.1.8", "EC3.2.1.37", "EC3.2.1.21",
                          "EC3.2.1.23", "EC3.2.1.55",
                          "EC5.3.1.12", "EC3.1.1.11", "EC3.2.1.82", "EC3.2.1.89", 
                          "EC3.2.1.172", "EC4.2.2.2"),
                     cat=c(rep("Cellulose degradation", 7),
                           rep("Cellulose or pectin degradation", 2),
                           rep("Pectin degradation", 6)))
ec_all$ec <- factor(ec_all$ec, levels=ec_all$ec)

# cereal fiber
# the 9 selected species:
# 1. species were associated with fiber intake and beta>0, to generate the pool_enzyme_ecs dataset
# 2. species were associated with cereal fiber or fruit fiber
# 3. ecs were pec_ec or cel_ec & species-specific ec were associated with cereal fiber or fruit fiber

ceraf0 <- ecsresults_INT %>% filter(exposure=="ceraf" & 
                                      species %in% c("Bifidobacterium_adolescentis",
                                                     "Bacteroides_eggerthii",
                                                     "Bacteroides_vulgatus",
                                                     "Bifidobacterium_pseudocatenulatum",
                                                     
                                                     "Dorea_longicatena",
                                                     "Eubacterium_eligens",
                                                     "Fusicatenibacter_saccharivorans",
                                                     "Gemmiger_formicilis",
                                                     "Ruminococcus_lactaris") &
                                      ec %in% ec_all$ec) %>% 
  mutate(ec=factor(ec, levels=ec_all$ec)) %>% arrange(ec)
ceraf0$enzy <- factor(ceraf0$enzy, levels=rev(unique(ceraf0$enzy)))
ceraf0$species <- gsub("_", " ", ceraf0$species)
ceraf0$species <- factor(ceraf0$species, 
                         levels=c("Bifidobacterium adolescentis",
                                  "Eubacterium eligens",
                                  "Bacteroides vulgatus",
                                  "Bacteroides eggerthii",
                                  "Fusicatenibacter saccharivorans",
                                  "Gemmiger formicilis",
                                  "Dorea longicatena",
                                  "Ruminococcus lactaris",
                                  "Bifidobacterium pseudocatenulatum"))
ceraf0 <- merge(ceraf0, ec_all, by="ec")  
ceraf0$star <- ifelse(ceraf0$qval.fdr<0.25, "*", "")

# fruit fiber
frtaf0 <- ecsresults_INT %>% filter(exposure=="frtaf" & 
                                      species %in% c("Bifidobacterium_adolescentis",
                                                     "Bacteroides_eggerthii",
                                                     "Bacteroides_vulgatus",
                                                     "Bifidobacterium_pseudocatenulatum",
                                                     
                                                     "Dorea_longicatena",
                                                     "Eubacterium_eligens",
                                                     "Fusicatenibacter_saccharivorans",
                                                     "Gemmiger_formicilis",
                                                     "Ruminococcus_lactaris") &
                                      ec %in% ec_all$ec) %>% 
  mutate(ec=factor(ec, levels=ec_all$ec)) %>% arrange(ec)
frtaf0$coef[frtaf0$coef>0.5] <- 0.5

frtaf0$enzy <- factor(frtaf0$enzy, levels=rev(unique(frtaf0$enzy)))
frtaf0$species <- gsub("_", " ", frtaf0$species)
frtaf0$species <- factor(frtaf0$species, 
                         levels=c("Bifidobacterium adolescentis",
                                  "Eubacterium eligens",
                                  "Bacteroides vulgatus",
                                  "Bacteroides eggerthii",
                                  "Fusicatenibacter saccharivorans",
                                  "Gemmiger formicilis",
                                  "Dorea longicatena",
                                  "Ruminococcus lactaris",
                                  "Bifidobacterium pseudocatenulatum"))
frtaf0 <- merge(frtaf0, ec_all, by="ec")  
frtaf0$star <- ifelse(frtaf0$qval.fdr<0.25, "*", "")

# make plot
catcol <- c("#d7191c77", "#ccae61", "#fdae6177")
names(catcol) <- c("Cellulose degradation", "Cellulose or pectin degradation", "Pectin degradation")
ec_all <- merge(ec_all, enzyme[!duplicated(enzyme$enzy),c("ec", "enzy")], by="ec", all.x = T)
ec_all$enzy[ec_all$ec=="EC2.4.1.281"] <-
  "EC2.4.1.281: 4-O-beta-D-mannosyl...phosphorylase"
ec_all$enzy[ec_all$ec=="EC3.2.1.55"] <-
  "EC3.2.1.55: Non-reducing end...arabinofuranosidase"
ec_all <- ec_all %>% arrange(ec)
ec_all$enzy <- factor(ec_all$enzy, levels=rev(ec_all$enzy))
ec_all$cat <- factor(ec_all$cat, levels=c("Cellulose degradation", 
                                          "Cellulose or pectin degradation", "Pectin degradation"))

p1_0 <- ggplot(ec_all, aes(x=1, y=enzy)) +
  geom_tile(aes(fill=cat)) +
  labs(fill="Enzyme category", x="") +
  scale_fill_manual(values=catcol) +
  scale_y_discrete(position = "right") +
  theme_bw() +
  theme(panel.border = element_rect(linewidth=1.5), 
        panel.grid = element_blank(),
        axis.text.y = element_text(color = "black", size = 12),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        axis.title =element_blank(),
        legend.position = "none")

p1_1 <- ggplot(aes(x=species, y=enzy), data=ceraf0) +
  geom_point(aes(fill=coef), shape=21, size=7) +
  labs(x="", y="") +
  scale_fill_gradient2(low="#08519c",
                       mid="white",
                       high="#bd0026",
                       midpoint = 0,
                       guide = "colourbar",
                       breaks=c(-0.5,0,0.5),
                       limits=c(-0.5,0.5)) +
  scale_x_discrete(position = "top") +
  geom_text(aes(label=star), color="black", size=7, nudge_y = -0.25) +
  theme_bw() +
  theme(panel.border = element_rect(linewidth=1.5), 
        panel.grid = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank(),
        axis.text.x = element_text(angle = 60, face="italic", colour = "black", size=12,
                                   hjust=0, vjust=0),
        legend.position = "none")

p2_1 <- ggplot(aes(x=species, y=enzy), data=frtaf0) +
  geom_point(aes(fill=coef), shape=21, size=7) +
  labs(x="", y="") +
  scale_fill_gradient2(low="#08519c",
                       mid="white",
                       high="#bd0026",
                       midpoint = 0,
                       guide = "colourbar",
                       breaks=c(-0.5,0,0.5),
                       limits=c(-0.5,0.5)) +
  scale_x_discrete(position = "top") +
  geom_text(aes(label=star), color="black", size=7, nudge_y = -0.25) +
  theme_bw() +
  theme(panel.border = element_rect(linewidth=1.5), 
        panel.grid = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank(),
        axis.text.x = element_blank(),
        legend.position = "none")

pdf("/n/home03/ywan/fiber/0430/fiber_ecs.pdf", 
    width = 7.1, height = 10.1, onefile = F) # Open a new pdf file
egg::ggarrange(p1_1, p1_0,
               p2_1, p1_0, 
               nrow = 2,
               widths = c(1,0.06))
dev.off() # Close the file

legend0 <- ggplot(ec_all, aes(x=1, y=enzy)) +
  geom_tile(aes(fill=cat)) +
  labs(fill="Enzyme category", x="") +
  scale_fill_manual(values=catcol) +
  theme_bw() +
  theme(panel.border = element_rect(linewidth=1.5), 
        panel.grid = element_blank(),
        legend.position = "right",
        legend.text = element_text(colour = "black", size = 12),
        legend.title = element_text(colour = "black", size = 12))
pdf("/n/home03/ywan/fiber/0430/fiber_ecs_lgd1.pdf", 
    width = 3, height = 2, onefile = F) # Open a new pdf file
grid::grid.draw(ggpubr::get_legend(legend0))
dev.off() # Close the file

legend1 <- ggplot(aes(x=species, y=enzy), data=ceraf0) +
  geom_point(aes(fill=coef), shape=21, size=7) +
  labs(x="", y="", fill="Fiber-enzyme\nassociation ") +
  scale_fill_gradient2(low="#08519c",
                       mid="white",
                       high="#bd0026",
                       midpoint = 0,
                       guide = "colourbar",
                       breaks=c(-0.5,0,0.5),
                       limits=c(-0.5,0.5)) +
  theme(legend.position = "top") +
  guides(fill = guide_colourbar(title.position="top",
                                title.hjust = 0.5,
                                title.theme = element_text(colour = "black", size = 16),
                                label.theme = element_text(colour = "black", size = 16),
                                barheight = 1.5,
                                barwidth = 10))
pdf("/n/home03/ywan/fiber/0430/fiber_ecs_lgd2.pdf", 
    width = 3, height = 2, onefile = F) # Open a new pdf file
grid::grid.draw(ggpubr::get_legend(legend1))
dev.off() # Close the file

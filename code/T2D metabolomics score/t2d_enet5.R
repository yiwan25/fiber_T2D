load("~/Fiber/metabolomics/datapool.RData")
all=pool2
gdata::keep(all, f, met275, sure=T)

library(glmnet)
library(survival)

set.seed(1234567)
ind <- sample(1:dim(all)[1], dim(all)[1])
aa <- split(ind, cut(seq_along(ind), 1000, labels = FALSE))

record_db2 = data.frame(ind=1:1000,
                        NoMetabs=NA,
                        Cindex=NA)
met_coef_db2 = data.frame(HMDB=colnames(all[,c(1:275)]))
score <- list()

for (i in 401:500) { 
  train=all[-aa[[i]], ]
  
  x=as.matrix(train[,c(1:275)])
  y=as.matrix(data.frame(time=train$tdb2, status=train$type2dbv))
  
  # Harell C index, a higher C index means better prediction performance
  Training_CV=cv.glmnet(x, y, family = "cox", alpha=0.5, type.measure = "C", nfolds = 10)
  lambda_1se_10F = Training_CV$lambda.1se
  
  tempcoef = as.vector(coef(Training_CV, s=lambda_1se_10F)[,])
  met_coef_db2[,dim(met_coef_db2)[2]+1]=tempcoef
  names(met_coef_db2)[dim(met_coef_db2)[2]] = paste("ind_",i,sep='')
  
  record_db2[i,"NoMetabs"] = sum(met_coef_db2[,dim(met_coef_db2)[2]]!=0)
  
  lambda_c=data.frame(lambda=Training_CV$lambda, c=Training_CV$cvm)
  record_db2[i,"Cindex"] = lambda_c[lambda_c$lambda==lambda_1se_10F,]$c
  
  test=all[aa[[i]], ]
  identical(colnames(test[, 1:275]), met_coef_db2$HMDB)
  
  tempscore = data.frame(matrix(nrow=length(aa[[i]]), ncol=2))
  colnames(tempscore) = c("newid", "score")
  tempscore$newid <- test$newid
  tempscore$score <- colSums(t(test[, 1:275])*tempcoef)
  score[[i]] <- tempscore
}

save.image(file="/udd/nhyiw/Fiber/metabolomics/score/t2d/t2d_enet5.RData")

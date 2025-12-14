#Code to run consensus clustering, sigclust, and limma

library("ConsensusClusterPlus")

d<-read.csv('../data/processed_data/m6a_for_ccp_top25stdmean_wbreast.csv', row.names="X")
d=data.matrix(d)
results = ConsensusClusterPlus(d, maxK=6, reps=1000, pItem=0.8, pFeature=1, title="ccp_m6a_wbreast", clusterAlg="hc", distance="pearson", seed=123, plot="png")

write.csv(results[[2]]$consensusClass, "../results/ccp_m6a_wbreast/wbreast_labels2.csv")
write.csv(results[[3]]$consensusClass, "../results/ccp_m6a_wbreast/wbreast_labels3.csv")
write.csv(results[[4]]$consensusClass, "../results/ccp_m6a_wbreast/wbreast_labels4.csv")
write.csv(results[[5]]$consensusClass, "../results/ccp_m6a_wbreast/wbreast_labels5.csv")
write.csv(results[[6]]$consensusClass, "../results/ccp_m6a_wbreast/wbreast_labels6.csv")

source("SigClust-v1.R")

k5<-read.csv("../data/processed_data/m6a_k5ccp_wbreast_for_sigclust.csv", row.names="X")
k5<-data.matrix(k5)
ps<-process_sigclust(k5, title = "m6a_sigclust_k5_wbreast")

library("limma")

m6a = read.csv("../data/processed_data/m6a_imputed_data_unfiltered.csv", row.names = "X")
design = read.csv("../data/processed_data/cluster_design_for_limma.csv", row.names = "X")


fit <- lmFit(m6a, design)

cont.matrix <-makeContrasts(X1vsAll = X1 - (X2+X3+X4+X5), levels = design)
fit2 <- contrasts.fit(fit, cont.matrix)
fit2 <- eBayes(fit2)


write.csv(topTable(fit2, adjust = "BH", p.value = 0.05, lfc = 1, coef = "X1vsAll", number = 20000), "../results/cluster_limma/1vsall.csv")


cont.matrix <-makeContrasts(X2vsAll = X2 - (X1+X3+X4+X5), levels = design)
fit2 <- contrasts.fit(fit, cont.matrix)
fit2 <- eBayes(fit2)


write.csv(topTable(fit2, adjust = "BH", p.value = 0.05, lfc = 1, coef = "X2vsAll", number = 20000), "../results/cluster_limma/2vsall.csv")


cont.matrix <-makeContrasts(X3vsAll = X3 - (X2+X1+X4+X5), levels = design)
fit2 <- contrasts.fit(fit, cont.matrix)
fit2 <- eBayes(fit2)


write.csv(topTable(fit2, adjust = "BH", p.value = 0.05, lfc = 1, coef = "X3vsAll", number = 20000), "../results/cluster_limma/3vsall.csv")


cont.matrix <-makeContrasts(X4vsAll = X4 - (X2+X3+X1+X5), levels = design)
fit2 <- contrasts.fit(fit, cont.matrix)
fit2 <- eBayes(fit2)


write.csv(topTable(fit2, adjust = "BH", p.value = 0.05, lfc = 1, coef = "X4vsAll", number = 20000), "../results/cluster_limma/4vsall.csv")

cont.matrix <-makeContrasts(X5vsAll = X5 - (X2+X3+X4+X1), levels = design)
fit2 <- contrasts.fit(fit, cont.matrix)
fit2 <- eBayes(fit2)


write.csv(topTable(fit2, adjust = "BH", p.value = 0.05, lfc = 1, coef = "X5vsAll", number = 20000), "../results/cluster_limma/5vsall.csv")


## Prepare TaqMan data

taq <- read.table("../data/TAQ.txt", stringsAsFactors=FALSE, header=TRUE, sep="\t")
## remove duplicates
taq <-  taq[!duplicated(taq[,1]),]
rownames(taq) <- taq[,1]
taq <- taq[,-1]

mean.a <- rowMeans(taq[, grepl("A", colnames(taq))])
mean.b <- rowMeans(taq[, grepl("B", colnames(taq))])

taq.dat <- cbind(mean.a, mean.b, logFC=log2(mean.b)- log2(mean.a))
save(taq.dat, file="../data/TaqManData.Rdata")

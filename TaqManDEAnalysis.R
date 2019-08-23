## Load DE results
load("../../ERCC_analysis/BGI_Lane01/.RData")

## Load TaqMan Data
load("../TaqManData.Rdata")
taq.data <- as.data.frame(taq.dat, stringsAsFactors=FALSE)


###########
## DESeq2
###########

deseq.taq <- merge(res.w.ercc1,
                   ##res.deseq$all.res, ## use this if using Results_tophat2.RData
                   taq.dat,
                   by.x=1, by.y='row.names')


## remove NaN
deseq.taq <- deseq.taq[-which(is.na(deseq.taq[,4])),]

## change infinite to +/- 1 from max or min
infinite.values <- deseq.taq[is.infinite(deseq.taq[,'log2FoldChange']), 'log2FoldChange'] 
infinite.values <-  sapply(infinite.values, function(x)
                           ifelse(sign(x)==1,  max(deseq.taq[is.finite(deseq.taq[,'log2FoldChange']),
                                        'log2FoldChange']) +1,
                                   min(deseq.taq[is.finite(deseq.taq[,'log2FoldChange']), 'log2FoldChange']) -1))

deseq.taq[is.infinite(deseq.taq[,'log2FoldChange']), 'log2FoldChange'] <- infinite.values



####################
## plot correlation
####################
PlotCorrelation <- function(dat, x.index, y.index, title, col, type='p'){
  ## Scatter plot of method logFC and TaqMan logFC
  ## type specifies correlation method
  ## p = Pearson (default), c= Concordance, r=RMS
  ## 
  require(MASS)
  require(epiR)
  datX <- dat[,x.index]
  datY <- dat[,y.index]

  ## Check type of correlation
  if(type == 'p'){ ## Pearson
    correlation <- cor(datX , datY)
  }else if(type == 'c'){ ## Concorance
    correlation <- epi.ccc(datX, datY)
  } else if (type == 'r'){ ## RMSD
    correlation <- sqrt(sum(((datX - datY)^2))/length(datX))
  } else{
    stop("Correlation type is unsepcified or unrecognized. must be 'p' for Pearson, 'c' Concordance or 'r' for RMSD")
  }

  
  plot(datX , datY,
       pch=19,
       col=col,
       main=title,
       sub=paste("n=",dim(dat)[1]),
       xlim=c(-15,15), ylim=c(-15,15),
       xlab=paste(title, 'logFC', sep=' '),
       ylab='TaqMan logFC')

  ## plot best fit line
  model <- rlm(datY ~ datX)
  abline(model, col='black')
  abline(0,1, col='red')

  ## add correlation values to title
  if(type == 'p'){
    mtext(paste("Pearson correlation", format(correlation,digits=3), sep=' '))
  }else if(type == 'c'){
    mtext(paste("Concordance correlation ", format(correlation$rho.c[,1], digits=3),
                " (95% CI ", format(correlation$rho.c[,2], digits=3), ' - ',
                format(correlation$rho.c[,3], digits=3), ")", sep=''))
  } else if (type == 'r'){
    mtext(paste("RMSD correlation ", format(correlation, digits=3), sep=''))
  }
  
  return(correlation)
}


#####################
#Lingqi script
#####################

pdf(paste("./TaqManAnalysis_", Sys.Date(), ".pdf", sep=''))

## plot Pearson correlations of uncentered data
p.cor.res <- PlotCorrelation(deseq.taq, 3, dim(deseq.taq)[2], "Deseq2 Pearson correlations","blue", type='p')
## Concorence correlation
c.cor.res <- PlotCorrelation(deseq.taq, 3, dim(deseq.taq)[2], "Deseq2 Concordence correlations","green", type='c')
## RMSD correlations
rms.cor.res <- PlotCorrelation(deseq.taq, 3, dim(deseq.taq)[2], "Deseq2 RMSD correlations","red", type='r')


names(rms.cor.res) <- "DESeq2"
## plot summary of correlations
summary.cor <- cbind(p.cor.res,
                     c.cor.res$rho.c[1],
                     #apply(c.cor.res, 2, function(x) x[["rho.c"]][[1]]),
                     rms.cor.res)
rownames(summary.cor) <- "DESeq2"
colnames(summary.cor) <- c("Pearson", "Concordance", "RMSD")

barplot(t(summary.cor), beside=T, legend=colnames(summary.cor),
        cex.names=0.75,
        ##ylim=c(0,1),
        main="Correlation summary")

barplot(rms.cor.res, col="blue", cex.name=.75,
        ylim=c(0, 2.5),
        ylab='RMSD from QRT-PCR log2 expression changes',
        main="RMSD correlation with TaqMan fold changes")

kLog2Cutoff <- 0.5
##############
## plot ROC
##############
PlotRocs <- function(dat, qval.index, logFC.index){
  require(pROC)
  outcome= rep(1, dim(dat)[1])
  outcome[abs(dat[,logFC.index]) <= kLog2Cutoff] =0
  
  roc <- plot.roc(outcome, dat[,qval.index],col="blue",
                    main="ROC of TaqMan data", ylim=c(0,1.05))
  mtext(paste("logFC cutoff= ", kLog2Cutoff, sep=''), side=3, padj=-1.75, cex=.8)

  return(roc)
}


roc.res.taq <- PlotRocs(deseq.taq, 7, dim(deseq.taq)[2])

legends <- paste("DESeq2", "AUC =", format(roc.res.taq$auc, digits=3), sep=' ')
legend("bottomright", legend=legends, col="blue", lwd=3, cex=.75, inset=c(0,0.03))

#########################
## Calculate AUCs
## by changing log2 cutoff
##########################
x_AUC <- function(dat, qval.index, logFC.index){
  ## calculate ROC
  ## return AUC vector for a range of logFC cutoffs
  require(pROC)
  auc.res <- matrix(nrow=length(seq(0.5,2,0.1)), ncol=1)
  
  ## logFC cutoff range
  cutoff <- seq(0.5,2,0.1)
  
  for(i in seq(1:length(cutoff))){
    outcome <- rep(1, dim(dat)[1])
    outcome[abs(dat[,logFC.index]) <= cutoff[i]] =0
    
    auc.res[i] <- roc(outcome, dat[,qval.index])$auc[[1]]
  }
  return(auc.res)
}

auc.res <- x_AUC(deseq.taq, 7, dim(deseq.taq)[2]) ## TaqMan logFC is last column


## plot AUCs
plot(seq(0.5,2,0.1), auc.res[,1], type='b', main="TaqMan AUCs",
     xlab="logFC cutoff values", ylab="AUC", col = "blue",
     ylim=c(0.8,1))

#for(i in seq(dim(auc.res)[2])){
#  lines(seq(0.5,2,0.1), auc.res[,i],
#        lwd=3, col=colr[i])
#}
legend("topleft", legend="DESeq2", col="blue", lwd=3,  cex=.75)


dev.off()
######################
#End of Lingqi script
######################




######################
## Center data
## by sample and gene
######################
CenterData <- function(dat){

  ## substract sample median values
  dat2=sweep(dat,2,apply(dat,2,median),"-")

  ## divide by sample sd.
  sd <- apply(as.matrix(dat2), 2, sd)
  dat3 <- dat2
  for(i in 1:ncol(dat2)){
    dat3[,i] <- dat2[,i]/sd[i]
  }

  ## substract gene median
  dat4=sweep(dat3,1,apply(dat3,1,median),"-")
  return(dat4)
}

###################
## Find common
## subset of genes
###################
x_intersect <- function(dat){
  ## recursive intersection function to
  ## find common element in list of elements
  if(length(dat) ==2){
    return (intersect(dat[[1]], dat[[2]]))
  }else{
    return (intersect(dat[[length(dat)-1]],x_intersect(dat[-(length(dat)-1)])))
  }
} 

FindCommonSubset <- function(dat){
  ## res <- intersect(dat[[1]],intersect(dat[[2]],intersect(dat[[3]],intersect(dat[[4]], dat[[5]]))))
  res <- x_intersect(dat)
  indicies <- sapply(dat, function(x) which(x %in% res))
  rownames(indicies) <- res
  return(indicies)
}

## list of DE methods logFC columns
## note that TaqMan logFC values are alway the last column
## in the matrix
logFC.index <- list(DESeq=4, edgeR=5,limmaQN=4,limmaVoom=4, PoissonSeq=5, CuffDiff=10, baySeq=2)

## extract the common genes
gene.names <- lapply(plot.dat, function(x) x[,1])
indicies <- FindCommonSubset(gene.names)

## generate a matrix of logFC
corr.dat <- sapply(seq(kNumOfMethods), function(i) plot.dat[[i]][indicies[,i],logFC.index[[i]]])
corr.dat <- cbind(corr.dat, taq.dat[rownames(indicies),"logFC"])
colnames(corr.dat) <- c(colnames(indicies), "TaqMan")

## Note: This centered data calculation
## is wrong. Centering should be done on
## normalized read counts not logFC values
## corr.dat <- CenterData(corr.dat)

## Color scheme by color pallette
## color.function <- colorRampPalette(c("orange", "gray", "red"), space='rgb')
## color.function <- colorRampPalette(c("orange", "gray", "#B51636"), space='rgb')
## colr <- color.function(kNumOfMethods)
colr <- c("#A6CEE3", "#1F78B4", "#B2DF8A", "#33A02C", "#FB9A99", "#E31A1C", "#FDBF6F")

plot2file=FALSE
if(plot2file){
  pdf(paste("../results/TaqManAnalysis_", Sys.Date(), ".pdf", sep=''))
}

## plot Pearson correlations of uncentered data
p.cor.res <- sapply(seq(kNumOfMethods), function(i) PlotCorrelation(plot.dat[[i]],
                                           logFC.index[[i]],
                                           dim(plot.dat[[i]])[2],
                                           names(plot.dat)[i],
                                           colr[i], type='p'))

## Concorence correlation
c.cor.res <- sapply(seq(kNumOfMethods), function(i) PlotCorrelation(plot.dat[[i]],
                                           logFC.index[[i]],
                                           dim(plot.dat[[i]])[2],
                                           names(plot.dat)[i],
                                           colr[i], type='c'))
## RMSD correlations
rms.cor.res <- sapply(seq(kNumOfMethods), function(i) PlotCorrelation(plot.dat[[i]],
                                           logFC.index[[i]],
                                           dim(plot.dat[[i]])[2],
                                           names(plot.dat)[i],
                                           colr[i], type='r'))

names(rms.cor.res) <- names(plot.dat)
## plot summary of correlations
summary.cor <- cbind(p.cor.res,
                     apply(c.cor.res, 2, function(x) x[["rho.c"]][[1]]),
                     rms.cor.res)
rownames(summary.cor) <- names(plot.dat)
colnames(summary.cor) <- c("Pearson", "Concordance", "RMSD")

barplot(t(summary.cor), beside=T, legend=colnames(summary.cor),
        cex.names=0.75,
        ##ylim=c(0,1),
        main="Correlation summary")

barplot(rms.cor.res, col=colr, cex.name=.75,
        ylim=c(0, 2.5),
        ylab='RMSD from QRT-PCR log2 expression changes',
        main="RMSD correlation with TaqMan fold changes")
        

###################################
## plot Sensitivity and Specificity
## at a fixed logFC cutoff and fixed adj.p.val
####################################
PlotQvalCorrelation <- function(dat, x.index, y.index, kLog2Cutoff, title, col){
  TP <- abs(dat[,x.index]) >= kLog2Cutoff & -log10(dat[,y.index]) >= 1.3
  
  ## sensitivity = TP/(TP+FN)
  sensitivity <- length(which(TP))/length(which(abs(dat[,x.index]) >= kLog2Cutoff))

  ## specificity =  TN/(TN + FP)
  TN <- abs(dat[,x.index]) < kLog2Cutoff & -log10(dat[,y.index]) < 1.3
  specificity <- length(which(TN))/length(which(abs(dat[,x.index]) < kLog2Cutoff))
  
  plot(dat[,x.index], -log10(dat[,y.index]),
       pch=19,
       col=ifelse(TP, col,ifelse(TN, 'orange','gray')),
       main=title,
       ## xaxt="n",
       xlab="TaqMan logFC",
       ylab="-10log(q-values)")
  abline(h=1.3, lty=2)
  abline(v=kLog2Cutoff, lty=3)
  abline(v=-kLog2Cutoff, lty=3)

  ## plot Sensitivity and sensitivity values
  ## and color-code the legend text to match the point colors
  positions <- par()$usr
  sensitivity <- format(sensitivity,digits=2)
  text(x=positions[1]*.95, y=positions[4]*.9,
       expression("Sensitivity= " * phantom(sensitivity)), adj=c(0,0))

  text(x=positions[1]*.95, y=positions[4]*.9,
       bquote(phantom("Sensitivity= ") * .(format(sensitivity,digits=2))), adj=c(0,0),
       col=c('blue'))

  text(x=positions[1]*.95, y=positions[4]*.85,
       expression("Specificity= " * phantom(specificity)), adj=c(0,0),
       col=c('black'))
  
  text(x=positions[1]*.95, y=positions[4]*.85,
       bquote(phantom("Specificity= ") * .(format(specificity,digits=2))), adj=c(0,0),
       col=c('orange'))
}

kLog2Cutoff <- 0.5 
## list of adj.pval columns or FDR
qval.index <- list(DESeq=3,  edgeR=4, limmaQN=3,limmaVoom=3, PoissonSeq=4, CuffDiff=13, baySeq=5)

## plot sensitivty and specificity
sapply(seq(kNumOfMethods), function(i) PlotQvalCorrelation(plot.dat[[i]],
                                               dim(plot.dat[[i]])[2],
                                               qval.index[[i]],
                                               kLog2Cutoff,
                                               names(plot.dat)[i],
                                               "blue"))


##############
## plot ROC
##############
PlotRocs <- function(i, dat, qval.index, logFC.index, color){
  require(pROC)
  outcome= rep(1, dim(dat)[1])
  outcome[abs(dat[,logFC.index]) <= kLog2Cutoff] =0
  if(i==1){
    roc <- plot.roc(outcome, dat[,qval.index],col=color,
                    main="ROC of TaqMan data", ylim=c(0,1.05))
    mtext(paste("logFC cutoff= ", kLog2Cutoff, sep=''), side=3, padj=-1.75, cex=.8)

  }else{
    roc <- lines.roc(outcome, dat[,qval.index], add=TRUE, col=color)
  }
  return(roc)
}


res <- lapply(seq(kNumOfMethods), function(i) PlotRocs(i, plot.dat[[i]],
                        qval.index[[i]],
                        dim(plot.dat[[i]])[2],
                        ##logFC.index[[i]],
                        colr[i]))

names(res) <- names(plot.dat)
legends <- lapply(seq(kNumOfMethods), function(i) paste(names(res)[i], "AUC =", format(res[[i]]$auc, digits=3), sep=' '))
legend("bottomright", legend=legends, col=colr, lwd=3, cex=.75, inset=c(0,0.03))

#########################
## Calculate AUCs
## by changing log2 cutoff
##########################
x_AUC <- function(i, dat, qval.index, logFC.index){
  ## calculate ROC
  ## return AUC vector for a range of logFC cutoffs
  require(pROC)
  auc.res <- matrix(nrow=length(seq(0.5,2,0.1)), ncol=1)

  ## logFC cutoff range
  cutoff <- seq(0.5,2,0.1)

  for(i in seq(1:length(cutoff))){
    outcome <- rep(1, dim(dat)[1])
    outcome[abs(dat[,logFC.index]) <= cutoff[i]] =0
   
    auc.res[i] <- roc(outcome, dat[,qval.index])$auc[[1]]
  }
  return(auc.res)
}

auc.res <- sapply(seq(kNumOfMethods), function(i) x_AUC(i, plot.dat[[i]],
                        qval.index[[i]],
                        dim(plot.dat[[i]])[2])) ## TaqMan logFC is last column

colnames(auc.res) <- names(plot.dat)

## plot AUCs
plot(seq(0.5,2,0.1), auc.res[,1], type='n', main="TaqMan AUCs",
     xlab="logFC cutoff values", ylab="AUC",
     ylim=c(0.8,1))

for(i in seq(dim(auc.res)[2])){
  lines(seq(0.5,2,0.1), auc.res[,i],
        lwd=3, col=colr[i])
}
legend("topleft", legend=colnames(auc.res), col=colr, lwd=3,  cex=.75)

if(plot2file){
  dev.off()
}


## load ERCC data
path="/Users/lingqiluo/Dropbox/UPMC/Secondary_analysis/DNAnexus/RNASeq_Benchmark_Dataset/Data\ Analysis/ERCC_analysis/BGI_Lane01/"
ercc <- read.table(paste(path,"ERCC_Controls_Analysis.txt",sep=""), stringsAsFactors=FALSE, header=T, sep="\t", row.names=2)

## Load DE results
load(paste(path,"deseq2/deseq2.RData",sep=""))

############
## DESeq2
###########
# Calculate the baseMean by group
baseMeanPerGroup <- sapply(levels(dds$group), function(lvl) rowMeans(counts(dds,normalized=TRUE)[,dds$group == lvl]))
colnames(baseMeanPerGroup) <- c("baseMeanA","baseMeanB")

#extract data.frame data from dds result
res.w.ercc<-as.data.frame(res)

#cbind a 'id' column
res.w.ercc <- cbind(id=rownames(res),res.w.ercc)

#merge ercc into res.w.ercc
#temp <- merge(res.w.ercc, ercc, by.x = 'row.names', by.y = 'row.names')

#merge baseMeanPerGroup into res.w.ercc
res.w.ercc1 <- merge(res.w.ercc, baseMeanPerGroup, by.x='row.names', by.y='row.names')
rownames(res.w.ercc1) <- res.w.ercc1[,'Row.names']
res.w.ercc1 <- res.w.ercc1[,!(names(res.w.ercc1) %in% c("Row.names"))]

# write 'res.w.ercc1' to a file
write.table(res.w.ercc1, file = paste(path,"DESeq2_outcome_BGI_Lane01.txt",sep = ""), sep = "\t", col.names = T, row.names = F)

#merge ercc into res.w.ercc1
deseq.ercc <- merge(res.w.ercc1, ercc, by.x = 'row.names', by.y = 'row.names')
rownames(deseq.ercc) <- deseq.ercc$id
deseq.ercc <- deseq.ercc[,!(names(deseq.ercc) == 'Row.names')]


# initiate a plot as pdf
pdf(paste("./ERCCAnalysis_", Sys.Date(), "_high_expression_gt_7attmol.ul.pdf", sep=''))

#pdf(paste("./ERCCAnalysis_", Sys.Date(), "_all_expression.pdf", sep=''))
####################
## plot correlation
####################

# build linear model to get both 'slope' and 'intercept'
BuildLinearModel <- function(dat, observed, expected){
  model <- lm(paste0("log2(", observed, ")", " ~ log2(", expected, ") "), dat )
  return(model)
}


# Get concordance correlation coefficient
BuildconcordanceModel <- function(dat, observed, expected){
  require(epiR)
  model <- epi.ccc(dat[,observed], dat[,expected])
  return(model)
}
#Define 'PlotCorrelation'
PlotCorrelation <- function(dat, x.index, y.index, title, col){
  # remove rows with 'NA' in the column 'log2FoldChange'
  dat <- dat[! is.na(dat$log2FoldChange), ]
  #print(dat) # works!
  # add 'minus sign' to 'log2FoldChange', because it is log ratio of B/A or Mix2/Mix1
  dat[,y.index] <- -dat[,y.index]
  #print(dat) works!
  print(dat[,x.index])
  print(dat[,y.index])
  
  plot(dat[,x.index], dat[,y.index],
  pch=19,
  col=col,
  main=title,
  xlab='Expected ERCC mix Log2 ratio',
  ylab=paste('Observed Log2 FC', sep=' '),
  xlim=c(-2.5, 4),
  ylim=c(-2.5, 4)
  )
  
  # Reference line in red
  abline(0,1,col="red")

  # Build model and extract parameters
  #model <- BuildLinearModel(dat, , expected)
  model <- lm(paste0(y.index," ~ ", x.index), dat)
  model_ccc <- BuildconcordanceModel(dat, x.index, y.index)
  
  intercept <- model$coefficients["(Intercept)"]
  slope <- model$coefficients[x.index]
  n <- summary(model)$df[2] + 1
  adj_R2 <- summary(model)$adj.r.squared
  
  # Add line and text
  abline(intercept, slope)
  # text position changed for all expression set
  text(-1.2,2.7,paste(paste0("Linear model parameters: \n n = ", n), paste0("R2 = ", round(adj_R2, 4)), paste0("Slope = ", round(slope, 2)), sep = "\n"))
  text(-1.2,3.8,paste("Concordance correlation ", format(model_ccc$rho.c[,1], digits=3),
                  " \n(95% CI ", format(model_ccc$rho.c[,2], digits=3), ' - ',
                  format(model_ccc$rho.c[,3], digits=3), ")", sep=''))
  # text position changed for high expression set (conc. > 7)
  #text(-2,2,paste(paste0("n = ", n), paste0("R2 = ", round(adj_R2, 4)), paste0("Slope = ", round(slope, 2)), sep = "\n"))
  
}

PlotCorrelationNormCount <- function(dat, sampleid){
  mix_idx = 2
  if (sampleid == "A") {
    mix_idx = 1
  }
  expected <- paste0("concentration.in.Mix.",mix_idx,"..attomoles.ul.")
  observed <- paste0("baseMean",sampleid)
  title <- paste0("DESeq2 ERCC Dose Response \n BGI Lane01 \n High Expression (Sample ",sampleid, ")")
  col <- '#A6CEE3'
  
  # Add 1 to columns 'baseMean*' and 'concentration.in.Mix.*..attomoles.ul.'
  dat[, expected] <- dat[,expected] + 1
  dat[, observed] <- dat[,observed] + 1
  
  # Generate plot
  plot(log2(dat[,expected]), log2(dat[,observed]),
       pch=19,
       col=col,
       main=title,
       xlab='ERCC mix Log2 concentration \n (attomoles/ul)',
       ylab=paste('Log2 Normalized Count', sep=' '),
       xlim=c(0, 20),
       ylim=c(0, 20)
       )
  
  # Build model and extract parameters
  model <- BuildLinearModel(dat, observed, expected)
  intercept <- model$coefficients["(Intercept)"]
  slope <- model$coefficients[paste0("log2(",expected,")")]
  n <- summary(model)$df[2] + 1
  adj_R2 <- summary(model)$adj.r.squared
  
  # Add line and text
  abline(intercept, slope)
  #abline(lm_SampleA$coefficients["(Intercept)"], lm_SampleA$coefficients["log2(concentration.in.Mix.1..attomoles.ul.)"])
  text(2,17.5,paste(paste0("n = ", n), paste0("R2 = ", round(adj_R2, 4)), paste0("Slope = ", round(slope, 2)), sep = "\n"))
  
}


# remove those spike-in RNAs with conc. less than 7 attomoles/ul (remove low expressed genes since general read depth of our dataset is low)
deseq.ercc.gt.7 <- deseq.ercc[deseq.ercc$concentration.in.Mix.1..attomoles.ul. > 7.0 & deseq.ercc$concentration.in.Mix.2..attomoles.ul. > 7.0, ]

#Make scatterplot for both log2 ratio and dose response
PlotCorrelation(deseq.ercc.gt.7, 'log2.Mix.1.Mix.2.', 'log2FoldChange', "DESeq2 ERCC Log2 \n BGI Lane01 \n High Expression (Sample A/Sample B)", '#A6CEE3')
PlotCorrelationNormCount(deseq.ercc.gt.7,"A")
PlotCorrelationNormCount(deseq.ercc.gt.7,"B")

#Make scatterplot for both log2 ratio and dose response
#PlotCorrelation(deseq.ercc, 'log2.Mix.1.Mix.2.', 'log2FoldChange', "DESeq2 ERCC Log2 \n BGI Lane01 \n (Sample A/Sample B)", '#A6CEE3')
#PlotCorrelationNormCount(deseq.ercc,"A")
#PlotCorrelationNormCount(deseq.ercc,"B")

dev.off()
#####################
#removing low count spike-in RNAs before plotting might help
#####################


####################
## plot ROCs
####################
require(pROC)

#########################
#Lingqi's code
#########################
PlotRocs <- function(dat, qval.index, logmix.index, title, color){
  outcome = rep(1,dim(dat)[1])
  outcome[dat[,logmix.index] == 0] = 0
  roc <- plot.roc(outcome, dat[,qval.index],col = color, main = title)
  roc <- lines.roc(outcome, dat[,qval.index], add=TRUE, col=color)
  return(roc)
}
pdf(paste("./ERCCAnalysis_", Sys.Date(), "_ROC.pdf", sep=''))
deseq.ercc.roc <- PlotRocs(deseq.ercc, 7, ncol(deseq.ercc), "ROC of ERCC spike-in data", '#A6CEE3')
legend("bottomright", legend=paste0("AUC = ",format(deseq.ercc.roc$auc, digits = 3)), col='#A6CEE3', lwd=2,  cex=.75)
dev.off()

pdf(paste("./ERCCAnalysis_", Sys.Date(), "_high_expression_ROC.pdf", sep=''))
deseq.ercc.gt.7.roc <- PlotRocs(deseq.ercc.gt.7, 7, ncol(deseq.ercc.gt.7),  "DESeq2 ROC of ERCC spike-in data\n(High expression)",'#A6CEE3')
legend("bottomright", legend=paste0("AUC = ",format(deseq.ercc.gt.7.roc$auc, digits = 3)), col='#A6CEE3', lwd=2,  cex=.75)

dev.off()

#rownames(res$de) <- res$de[,'id']
#res.deseq$de <- res.deseq$de[,c('id', 'pval', 'padj', 'log2FoldChange', 'baseMeanA', 'baseMeanB')]

#deseq.taq <- merge(res.deseq$de,
                   ## res.deseq$all.res,
 #                  ercc,
  #                 by.x='row.names', by.y='row.names')

## change infinite to +/- 1 from max or min
#infinite.values <- deseq.taq[is.infinite(deseq.taq[,'log2FoldChange']), 'log2FoldChange'] 
#infinite.values <-  sapply(infinite.values, function(x)
 #                          ifelse(sign(x)==1,  max(deseq.taq[is.finite(deseq.taq[,'log2FoldChange']),
                                  #      'log2FoldChange']) +1,
  #                                 min(deseq.taq[is.finite(deseq.taq[,'log2FoldChange']), 'log2FoldChange']) -1))

#deseq.taq[is.infinite(deseq.taq[,'log2FoldChange']), 'log2FoldChange'] <- infinite.values
#plot.dat["DESeq"] <- list(deseq.taq)





if(FALSE){
## Boxplot of log2.mix1mix2=0 q-values 
## collect zero p-values for log2.Mix.1.Mix.2.
qval.index <- list(DESeq=4, edgeR=4, limmaQN=3,limmaVoom=3, PoissonSeq=4,
                   CuffDiff=13, baySeq=5, NOISeq=6)
qval.0.mix.dat <- lapply(seq(kNumOfMethods), function(i) plot.dat[[i]][
                                                 plot.dat[[i]][,dim(plot.dat[[i]])[2]] == 0,
                                                 qval.index[[i]]
                                                 ])


names(qval.0.mix.dat) <- names(plot.dat)
## boxplot(qval.0.mix.dat, ylab='Adjusted p-values', col='bisque',
##         main='ERCC Non differentiated spike-ins adjusted p-values')

require(beeswarm)
myCol <- lapply(seq(kNumOfMethods), function(i) ifelse(qval.0.mix.dat[[i]] <=0.05, 'gray', colr[i]))
names(myCol) <- names(plot.dat)
beeswarm(qval.0.mix.dat, pch=19, add=FALSE, pwcol=myCol,
         main="ERCC Non differentiated spike-ins adjusted p-values")
abline(h=0.05, lwd=2, lty=3)


## get the fraction of false positive q-values
print(sapply(qval.0.mix.dat, function(x) table(x<0.05)))
}


library(GEOquery)
library(wateRmelon)
library(IlluminaHumanMethylation450k.db)
library(lattice)
library(reshape2)
library(limma)
library(gplots)

if (file.exists("methyl_ALL.Rdata")) {
  # if previously downloaded
  load("methyl_ALL.Rdata")
} else {
  # if downloading for the first time
  GSE39141 <- getGEO("GSE39141")
  show(GSE39141)  ## 33 samples (29 ALL and 4 healthy B cells)
  GSE42865 <- getGEO("GSE42865")  # took ~2 mins for JB
  show(GSE42865)  ## 16 samples (9 healthy cells B cells and 7 other cells)
  
  # Extract expression matrices (turn into data frames at once)
  ALL.dat <- as.data.frame(exprs(GSE39141[[1]]))
  CTRL.dat <- as.data.frame(exprs(GSE42865[[1]]))
  
  # Obtain the meta-data for the samples and rename them perhaps?
  ALL.meta <- pData(phenoData(GSE39141[[1]]))
  CTRL.meta <- pData(phenoData(GSE42865[[1]]))
  
  # create some labels
  ALL.meta$Group <- c(rep("ALL", 29), rep("HBC", 4))
  ## ALL: Case; HBC: Healthy B Cells
  
  # Subset both meta-data and data for control (healthy) donors
  CTRL.meta <- droplevels(subset(CTRL.meta, grepl("Healthy donor", characteristics_ch1.1)))
  CTRL.dat <- subset(CTRL.dat, select = as.character(CTRL.meta$geo_accession))
  
  # Rename variables
  names(ALL.dat) <- paste(ALL.meta$Group, gsub("GSM", "", names(ALL.dat)), 
                          sep = "_")
  names(CTRL.dat) <- paste("HBC", gsub("GSM", "", names(CTRL.dat)), sep = "_")
  
  # save the data to avoid future re-downloading
  save(ALL.dat, CTRL.dat, ALL.meta, CTRL.meta, file = "methyl_ALL.Rdata")
}

#  distribution of Beta
meandat <- data.frame(Beta = c(rowMeans(ALL.dat,na.rm=T), rowMeans(CTRL.dat,na.rm=T)), case = c(rep("ALL", length(rowMeans(ALL.dat,na.rm=T))), rep("CTRL", length(rowMeans(CTRL.dat,na.rm=T)))))
densityplot(~Beta,meandat, groups = case, grid = TRUE, plot.points = FALSE, auto.key = TRUE,main="Average Beta value density of two experiments")
# distribution of Beta values are not the same between different experiments. 
# Normalization
# combine data from two experiments into one matrix, each column represents
# beta values of one sample
beta.matrix <- as.matrix(cbind(ALL.dat, CTRL.dat))
str(beta.matrix, max.level = 0)
# quantile normalization
system.time(beta.norm <- betaqn(beta.matrix))
# Density of beta values before and after nomalizaion
dat.probeMeans <- c(rowMeans(beta.norm[, 1:ncol(ALL.dat)], na.rm = TRUE),rowMeans(beta.norm[, ncol(ALL.dat):ncol(beta.norm)], na.rm = TRUE)) 
plotNorm <-rbind(data.frame(meandat, Norm = "Before"),data.frame(Beta = dat.probeMeans, case = rep(c('ALL', 'CTRL'), each = nrow(ALL.dat)),Norm = "After"))
plotNorm$Norm <- factor(plotNorm$Norm, levels = c("Before", "After"))
densityplot(~Beta | Norm,groups = case,data=plotNorm,grid = TRUE, plot.points = FALSE, auto.key = TRUE,main="Average Beta value before and after normalization")
# M values
M.norm <- beta2m(beta.norm)
# CpG Islands
# Extracting probe ID to CpG islands association
cginame <- as.data.frame(IlluminaHumanMethylation450kCPGINAME)
names(cginame) <- c("Probe_ID", "cginame")
rownames(cginame) <- cginame$Probe_ID
length(levels(factor(cginame$cginame)))  # No. of CGIs
# restrict probes to those within CGIs
beta.inCGI <- beta.norm[cginame$Probe_ID, ]
M.inCGI <- M.norm[cginame$Probe_ID, ]
nrow(M.inCGI)  # No. of probes within CGIs
# aggregate probes to CGIs
beta.CGI <- aggregate(beta.inCGI, by = list(cginame$cginame), mean, na.rm = T)
rownames(beta.CGI) <- beta.CGI[, "Group.1"]
beta.CGI <- subset(beta.CGI, select = -Group.1)
str(beta.CGI, max.level = 0)
M.CGI <- aggregate(M.inCGI, by = list(cginame$cginame), mean, na.rm = T)
rownames(M.CGI) <- M.CGI[, "Group.1"]
M.CGI <- subset(M.CGI, select = -Group.1)
str(M.CGI, max.level = 0)
M.CGI.tall <- melt(t(M.CGI), value.name = 'M', varnames = c('Sample', 'CGI'))
M.CGI.tall$Group <- gsub("_[0-9]+", "", M.CGI.tall$Sample)
head(M.CGI.tall)
colr <- rep("red",nrow(M.CGI.tall))
colr[which(M.CGI.tall$Group == "HBC")] <- "green"
suppressWarnings(boxplot(M~Sample,data=M.CGI.tall,ylab="M-values", xlab="Samples", main="Distribution of CGI M values",col=colr))
legend("bottomright", c("ALL", "HBC"), col = c("red", "green"), pch = 15)
# Red represents ALL, green represents HBC

# Differential methylation analysis with limma
design <- data.frame(Group = relevel(factor(gsub("_[0-9]+", "", colnames(M.CGI))), ref = "HBC"), row.names = colnames(M.CGI))
str(design)
DesMat <- model.matrix(~Group, design)
DMRfit <- lmFit(M.CGI, DesMat)
DMRfitEb <- eBayes(DMRfit)
cutoff <- 0.01
DMR <- topTable(DMRfitEb, coef = "GroupALL", number = Inf, p.value = cutoff)
head(DMR)  # top hits 
nrow(DMR)
# using a cutoff of FDR = 0.01, we identified 4115 CGIs 

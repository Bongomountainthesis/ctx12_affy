###dont know what we're doing but following instructions from 

#http://bioinf.wehi.edu.au/marray/ibc2004/lab2/lab2.html
#http://www.bioconductor.org/help/course-materials/2006/biocintro_april/thurs/affy/Affy.pdf

stringsAsFactors=FALSE

library(affy)
library(limma)
library(biomaRt)

library(affydata)

files <- list.celfiles(path = "/space/angela/ctx12_affy/data", full.names=TRUE)

data <- ReadAffy(filenames = files)

#do some QC stuff

postscript(file = "results/matplot_prenorm.ps",horizontal=FALSE)
matplot(pm(data),type="l", xlab ="Probe No.", ylab ="PM Probe intensity")
dev.off()

postscript(file = "results/matplot_t_prenorm.ps",horizontal=FALSE)
matplot(t(pm(data)),type="l", xlab ="Probe No.", ylab ="PM Probe intensity")
dev.off()

postscript(file = "results/histogram_prenorm.ps",horizontal=FALSE)
hist(data)
dev.off()

postscript(file = "results/boxplot_prenorm.ps",horizontal=FALSE)
boxplot(data)
dev.off()

#check phenotype data

pData(data)

#background correct data

#data.bg.rma <- bg.correct(data,method = "rma")

#normalise data

#data.norm <- normalize(data, method = "quantiles")

#use expresso??

data.expresso <- expresso(data, 
			normalize.method = "quantiles", 
			bgcorrect.method = "rma",
			pmcorrect.method = "pmonly",
			summary.method = "avgdiff"
			)

#get expression matrix and do some stuff

E <- exprs(data.expresso)

design<-matrix(0,nrow=(ncol(E)), ncol=3)
colnames(design) <- c("BMP4","CTX12","FBS")
rownames(design) <- colnames(E)
design[1:3,1] <- 1
design[4:6,2] <- 1
design[7:9,3] <- 1

cont.matrix<-makeContrasts(BMP4vsCTX12=BMP4-CTX12,
			   FBSvsCTX12=FBS-CTX12,
                           BMP4vsFBS=BMP4-FBS,
                           levels=design)

fit<-lmFit(E, design)
fit<-contrasts.fit(fit, cont.matrix)
ebFit<-eBayes(fit)

library(mouse4302.db)

ids = rownames(E)
symbol <- mget(ids, mouse4302SYMBOL, ifnotfound = NA)
sum(sapply(symbol, length) > 1)
symbol <- as.character(symbol)
length(which(symbol=="NA"))

ensembl = mget(ids, mouse4302ENSEMBL, ifnotfound = NA)
length(which(is.na(ensembl)))

sum(sapply(ensembl, length) > 1)
crosshyb <- which(( sapply(ensembl, length) ) > 1) 
length(crosshyb)
ensembl[crosshyb] <- NA 
ensembl <- as.character(ensembl)
ensembl[ensembl=="NA"] = NA
length(which(is.na(ensembl)))

ensmart <- useMart("ensembl", dataset="mmusculus_gene_ensembl")























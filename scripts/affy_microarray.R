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

#save stuff
save(data.expresso, file = "results/data_quantile.RData")

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

filters <- "ensembl_gene_id"
values <- ensembl[!is.na(ensembl)]
attributes <- c("ensembl_gene_id", "chromosome_name","start_position", "end_position", "strand", "description")
ens.anno <- getBM(filters=filters, values=values, attributes=attributes, mart=ensmart)
rownames(ens.anno)<-ens.anno[,1]

anno <- data.frame(
                   ID = as.character(ids),
                   EnsemblID = ensembl,
                   symbol=symbol,
                   ens.anno[ensembl,],
                   stringsAsFactors=F
              )
rownames(anno) <- anno[,"ID"]

ebFit$genes = anno 

#BMP4vsCTX12=BMP4-CTX12
#FBSvsCTX12=FBS-CTX12
#BMP4vsFBS=BMP4-FBS

f.test<- topTable(ebFit, number=nrow(E))
rownames(f.test)<-f.test$ID
f.test<-f.test[order(f.test[,"P.Value"]),]
write.csv(f.test,"results/f_test.csv",row.names=F)

BMP4vsCTX12 <- topTable(ebFit, coef=1, adjust="BH", number=nrow(E))
rownames(BMP4vsCTX12) <- BMP4vsCTX12$ID
BMP4vsCTX12 <- BMP4vsCTX12[order(abs(BMP4vsCTX12[,"logFC"]), decreasing=TRUE),]
BMP4vsCTX12 <- BMP4vsCTX12[!duplicated(BMP4vsCTX12[,"symbol"]),]
BMP4vsCTX12 <- BMP4vsCTX12[order(BMP4vsCTX12[,"adj.P.Val"],decreasing=FALSE),]
write.csv(BMP4vsCTX12,"results/BMP4vsCTX12.csv",row.names=F)

FBSvsCTX12<-topTable(ebFit, coef=2, adjust="BH", number=nrow(E))
rownames(FBSvsCTX12)<-FBSvsCTX12$ID
FBSvsCTX12 <- FBSvsCTX12[order(abs(FBSvsCTX12[,"logFC"]), decreasing=TRUE),]
FBSvsCTX12 <- FBSvsCTX12[!duplicated(FBSvsCTX12[,"symbol"]),]
FBSvsCTX12 <- FBSvsCTX12[order(FBSvsCTX12[,"adj.P.Val"],decreasing=FALSE),]
write.csv(FBSvsCTX12,"results/FBSvsCTX12.csv",row.names=F)

BMP4vsFBS <-topTable(ebFit, coef=3, adjust="BH", number=nrow(E))
rownames(BMP4vsFBS)<-BMP4vsFBS$ID
BMP4vsFBS <- BMP4vsFBS[order(abs(BMP4vsFBS[,"logFC"]), decreasing=TRUE),]
BMP4vsFBS <- BMP4vsFBS[!duplicated(BMP4vsFBS[,"symbol"]),]
BMP4vsFBS <- BMP4vsFBS[order(BMP4vsFBS[,"adj.P.Val"],decreasing=FALSE),]
write.csv(BMP4vsFBS,"results/BMP4vsFBS.csv",row.names=F)























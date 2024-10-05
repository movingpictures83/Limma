library("derfinder")
library("derfinderData")
library("GenomicRanges")
library("knitr")



#####################################################################################
# LIMMA
input <- function(inputfile) {
regionMat <<- readRDS(inputfile)
}

run <- function() {}

output <- function(outputfile) {
pheno <- subset(brainspanPheno, structure_acronym == "AMY")
p <- pheno[, -which(colnames(pheno) %in% c(
    "structure_acronym",
    "structure_name", "file"
))]
mod <- model.matrix(~ pheno$group + pheno$gender)
mod0 <- model.matrix(~ pheno$gender)
transformedCov <- log2(regionMat$chr21$coverageMatrix + 32)
library("limma")
fit <- lmFit(transformedCov, mod)
fit0 <- lmFit(transformedCov, mod0)
getF <- function(fit, fit0, theData) {
    rss1 <- rowSums((fitted(fit) - theData)^2)
    df1 <- ncol(fit$coef)
    rss0 <- rowSums((fitted(fit0) - theData)^2)
    df0 <- ncol(fit0$coef)
    fstat <- ((rss0 - rss1) / (df1 - df0)) / (rss1 / (ncol(theData) - df1))
    f_pval <- pf(fstat, df1 - df0, ncol(theData) - df1, lower.tail = FALSE)
    fout <- cbind(fstat, df1 - 1, ncol(theData) - df1, f_pval)
    colnames(fout)[2:3] <- c("df1", "df0")
    fout <- data.frame(fout)
    return(fout)
}
ff <- getF(fit, fit0, transformedCov)
limma <- regionMat$chr21$regions
limma$fstat <- ff$fstat
limma$pvalue <- ff$f_pval
limma$padj <- p.adjust(ff$f_pval, "BH")
#limma
write(table(limma$padj < 0.05), outputfile)
}
#####################################################################################


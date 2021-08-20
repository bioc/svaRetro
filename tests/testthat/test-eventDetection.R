context('event detection functions')
colo829 <- readVcf(system.file("extdata", "diploidSV.vcf", package = "svaRetro"))

#RT detection
#vcf <- VariantAnnotation::readVcf(vcf.file, "hg19")
#gr <- breakpointRanges(vcf, nominalPosition=TRUE)
#rt <- rtDetect(gr, genes, maxgap=30, minscore=0.6)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
genes <- TxDb.Hsapiens.UCSC.hg19.knownGene

test_that("RT detection returns two Granges", {
    gr <- StructuralVariantAnnotation::breakpointRanges(colo829, nominalPosition=TRUE)
    rt <- rtDetect(gr, genes, maxgap=50, minscore=0.3)
    if (length(rt)>0) {
        expect_equal(rep(2, length(rt)), unname(sapply(rt, length)))
    }
})

# this one returns "Error: txs should be a list object" on 
# .txs2genesym(rt.gr$txs) at eventDetection.R#135
# 21/AUG update: this error has been fixed

test_that("RT detection returns two Granges", {
    gr <- StructuralVariantAnnotation::breakpointRanges(colo829, nominalPosition=TRUE)
    rt <- rtDetect(gr, genes, maxgap=50, minscore=0.8)
    if (length(rt)>0) {
        expect_equal(rep(2, length(rt)), unname(sapply(rt, length)))
    }
})


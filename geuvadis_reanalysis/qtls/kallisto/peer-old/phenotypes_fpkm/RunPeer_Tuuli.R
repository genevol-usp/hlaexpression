##run PEER


arguments <- commandArgs(trailingOnly = TRUE)

inputfile <- arguments[1]
covariatefile <- arguments[2] #can be "-"
number_of_covs <- as.numeric(arguments[3])
number_of_iterations <- as.numeric(arguments[4])
output_prefix <- arguments[5]

#inputfile <- "YRI89.ExonQuantCount.45N.50F.samplename.chr20.tab"
#covariatefile <- "YRI89_mRNA_SeqLabNumber.txt"
#number_of_covs <- 20
#number_of_iterations <- 100
#output_prefix <- "YRI89.ExonQuantCount.45N.50Ftest"

library(peer)
expr = read.table(inputfile, header=FALSE)
model = PEER()

if (covariatefile != "-") {
        covs = read.table(covariatefile,header=FALSE) #NxC matrix
        PEER_setCovariates(model, as.matrix(covs))
}
PEER_setPhenoMean(model,as.matrix(expr))
PEER_setAdd_mean(model, TRUE)
PEER_setNk(model,number_of_covs)
PEER_setNmax_iterations(model, number_of_iterations)

PEER_update(model)

factors = PEER_getX(model)
weights = PEER_getW(model)
precision = PEER_getAlpha(model)
residuals = PEER_getResiduals(model)

write.table(factors, paste(output_prefix, "_factors.txt", sep=""), quote=F, sep="\t")
write.table(weights, paste(output_prefix, "_weights.txt", sep=""), quote=F, sep="\t")
write.table(precision, paste(output_prefix, "_precision.txt", sep=""), quote=F, sep="\t")
write.table(residuals, paste(output_prefix, "_residuals.txt", sep=""), quote=F, sep="\t")

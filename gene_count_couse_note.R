brary(edgeR)
# Read in gene count data matrix, set group info, filter expr and estimate dispersion.
df <- read.csv(file="example_gene_count_matrix.csv", row.names=1)
df[df==0] <- NA
df <- df[complete.cases(df), ]

# Set group info
group <- factor(rep(c('C','T'), times=c(3,3)))
obj <- DGEList(counts=df, group=group)
obj$samples              #check sample info 

# Filter low expr tags
filter <- filterByExpr(obj)
obj <- obj[filter, , keep.lib.sizes=F]

# Estimate data edispersion
design <- model.matrix(~0+group, data=obj$samples)
colnames(design) <- levels(obj$samples$group)
obj <- estimateDisp(obj,design)
obj$common.dispersion
plotBCV(obj)
plotMDS(obj, top=500, method = "bcv", col=c(rep(1:2, each=3)))

# To perform clasic exact test
et <- exactTest(obj, pair=c("C","T"))
topTags(et)
plotMD(et)
abline(h=c(-1,1), col="blue")
summary(decideTests(et))

#To perform quasi-likelihood F-tests (taking acount for technical and biological variations):
fit <- glmQLFit(obj, design, robust=T)
plotQLDisp(fit)
qlf <- glmQLFTest(fit,coef=2)
togTags(qlf)
summary(decideTests(qlf))
plotMD(qlf)
abline(h=c(-1,1), col="blue")
with(qlf$table, plot(logCPM,logFC,pch=16,cex=0.3))

#To perform likelihood ratio tests:
fit <- glmFit(obj,design)
lrt <- glmLRT(fit,coef=2)
topTags(lrt)


## When many more genes are in the background than in the gene set, p values of Wilcoxon test are lower than that of t-test.
getP <- function(seed, bgSize=10) {
    set.seed(seed)
    s1 <- rnorm(bgSize)
    s2 <- rnorm(10, mean=2)
    
    tP <- t.test(s2, s1, alternative="greater")$p.value
    tW <- wilcox.test(s2, s1, alternative="greater")$p.value
    return(c(tP, tW))
}

testTests <- function(N=1000, bgSize=10) {
    ps <- t(sapply(1:N, getP, bgSize=bgSize))
    colnames(ps) <- c("t-test", "wilcoxon")
    return(ps)
}

eqSize <- testTests(bgSize=10)
ieqSize <- testTests(bgSize=1000)
cols <- c("lightblue", "orange")

par(mfrow=c(1,2), mar=c(3,3,3,1), mgp=c(2,1,0))
boxplot(log10(eqSize), col=cols, ylab="log10P", main="Similar size")
boxplot(log10(ieqSize), col=cols, ylab="log10P", main="Bg>>gene set")
dev.print(pdf, "t-wilcoxon-size.pdf", width=6, height=4)

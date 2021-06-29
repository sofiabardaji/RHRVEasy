library("ggplot2")
library("dlstats")

x <- cran_stats(c( "RHRV"))

if (!is.null(x)) {
  head(x)
  ggplot(x, aes(end, downloads, group=package, color=package)) +
    geom_line() + geom_point(aes(shape=package))
}

x <- cran_stats(c( "RHRV"),use_cache = FALSE)
x

for(i in 1:39){
install.packages("RHRV")
}

a=RHRVEasy(folders =c("C:\\rrs\\normal3","C:\\rrs\\anormal33"), correction = TRUE, method = "none", nonLinear = TRUE)


a1=RHRVEasy(folders =c("C:\\rrsc\\normal","C:\\rrsc\\chf"), nonLinear = TRUE)
a2=RHRVEasy(folders =c("C:\\rrsc\\normal","C:\\rrsc\\chf"), correction = TRUE, method = "bonferroni")
a3=RHRVEasy(folders =c("C:\\rrsc\\normal","C:\\rrsc\\chf"), correction = TRUE, method = "holm")
a4=RHRVEasy(folders =c("C:\\rrsc\\normal","C:\\rrsc\\chf"), correction = TRUE, method = "hochberg")
#a5=RHRVEasy(folders =c("C:\\rrsc\\normal","C:\\rrsc\\chf"), correction = TRUE, method = "hommel")
a6=RHRVEasy(folders =c("C:\\rrsc\\normal","C:\\rrsc\\chf"), correction = TRUE, method = "fdr")
a7=RHRVEasy(folders =c("C:\\rrsc\\normal","C:\\rrsc\\chf"), correction = TRUE, method = "BY")

ac=RHRVEasy(folders =c("C:\\rrsc\\normal","C:\\rrsc\\chf"))
a2c=RHRVEasy(folders =c("C:\\rrsc\\normal","C:\\rrsc\\chf"), correction = TRUE)
a3c=RHRVEasy(folders =c("C:\\rrsc\\normal","C:\\rrsc\\chf"), correction = TRUE, method = "fdr")

cbind(a1$pValues,a2$pValues,a3$pValues,a4$pValues,a6$pValues,a7$pValues)

mapply('-', a1$pValues, a2$pValues, SIMPLIFY = FALSE)

cbind(a$pValues,ac$pValues)
cbind(a2$pValues,a3c$pValues)
cbind(a3$pValues,a3c$pValues)

a2=RHRVEasy(folders =c("C:\\rrs\\normal","C:\\rrs\\chf","C:\\rrs\\normal2"))
a21=RHRVEasy(folders =c("C:\\rrs\\normal","C:\\rrs\\chf","C:\\rrs\\normal2"), correction = TRUE)
a22=RHRVEasy(folders =c("C:\\rrs\\normal","C:\\rrs\\chf","C:\\rrs\\normal2"), correction = TRUE, method = "fdr")

b=RHRVEasy(folders =c("C:\\rrs\\normal","C:\\rrs\\chf"), correction = TRUE)
c=RHRVEasy(folders =c("C:\\rrs\\normal","C:\\rrs\\chf"), correction = TRUE, method = "fdr")

a2=RHRVEasy(directorios =c("C:\\rrs\\normal","C:\\rrs\\chf"),  correction = TRUE)
b2=RHRVEasy(control="C:\\rrs\\normal",case="C:\\rrs\\chf", correction = FALSE)
c2=RHRVEasy(control="C:\\rrs\\normal",case="C:\\rrs\\chf", correction = TRUE, method = "fdr")

a=read.csv("a.csv", sep=";")
write.csv(x,file = "dow.csv")















setwd("Z:\\articulos\\EPNEAC\\R")
x2 = read.csv("2.csv")
apply(x2, MARGIN=2, FUN=shapiro.test)

wilcox.test(x=x2$X1, y = x2$X2, exact = F, alternative ="greater")
wilcox.test(x=x2$X2, y = x2$X4, exact = F, alternative ="greater")
wilcox.test(x=x2$X4, y = x2$X6, exact = F, alternative ="greater")
wilcox.test(x=x2$X6, y = x2$X8, exact = F, alternative ="greater")
wilcox.test(x=x2$X8, y = x2$X10, exact = F, alternative ="greater")
wilcox.test(x=x2$X10, y = x2$X12, exact = F, alternative ="greater")

friedman.test(as.matrix(x2))
 t =posthoc.kruskal.nemenyi.test(x2)


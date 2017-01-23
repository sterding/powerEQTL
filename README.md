# powerEQTL
R package for power analysis for eQTL study

== TOADD == 

* Below is example code using the package to generate the similar power analysis curve as GTEx study[^1]

[^1]: Lonsdale J and Thomas J, et al. The Genotype-Tissue Expression (GTEx) project. Nature Genetics, 45:580-585, 2013.


```{r}
library('powerEQTL')
# sample size
N <- c(50,100,150,200,250,300)
nn <- length(N)

# MAF
MAF <- seq(0.5,20,0.1)/100
nq <- length(MAF)

# number of SNPs tested (10 SNPs per gene x 24 genes in total)
nSNP=240

# significant level (FP)
a=0.05

# obtain power
power_unbalanced <- array(numeric(nn*nq), dim=c(nn,nq))
for (i in 1:nn){
    for (j in 1:nq){
        # unbalanced
        result <- powerEQTL.ANOVA(MAF=MAF[j], 
                                  typeI=a, 
                                  nTests=nSNP, 
                                  myntotal=N[i], 
                                  mystddev=0.13,
                                  deltaVec = c(0.13, 0.13),
                                  verbose = F)
        power_unbalanced[i,j] <-result;
    }
}

# set up graph
xrange <- range(MAF*100)
yrange <- c(0:1)
colors <- rainbow(length(N))
plot(xrange, yrange, log='x', type="n",
     xlab="MAF (%)",
     ylab="Power",
     main="Power Estimation for eQTL Studies\nSig=0.05, nSNP=240 (one-way unbalanced ANOVA)")

abline(v=0, h=seq(0,1,.1), lty=2, col="grey89")
abline(h=0, v=c(1:10), lty=2,col="grey89")

# add power curves
for (i in 1:nn){
    lines(MAF*100, power_unbalanced[i,], type="l", lwd=4, col=colors[i])
}

legend("topleft", title="Sample size (n)", as.character(N),fill=colors, bty='n')
```
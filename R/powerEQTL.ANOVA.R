# created on Dec. 11, 2016
# (1) allow mu2-mu1 not equal to mu3-mu2,
#   where mu1 is the mean expression level for mutation homozygote,
# mu2 is the mean expression level for heterozygote,
#  and mu3 is the mean expression level for wildtype homozygote
#
# created on Dec. 8, 2016
#  (1) power calculation for eQTL analysis based on ANOVA or simple linear regression
#

# MAF - minor allele frequency
# typeI - type I error rate
# nTests - number of tests
# myntotal - total number of subjects
# mystddev - standard deviation of gene expression levels
#            (assume each group of subjects has the same mystddev)
# deltaVec - mean difference of gene expression levels between groups
#            deltaVec[1]=mu2-mu1 and deltaVec[2]=mu3-mu2
# verbose - flag indicating if we should output intermediate results

powerEQTL.ANOVA=function(MAF,
                   typeI=0.05,
                   nTests=200000,
                   myntotal=200,
                   mystddev=0.13,
                   deltaVec=c(0.13, 0.13),
                   verbose=TRUE)
{
  if(length(deltaVec)!=2)
  {
    stop("'deltaVec' has 2 and only 2 elements")
  }
  
  gm1 = -deltaVec[1]
  gm2 = 0
  gm3 = deltaVec[2]

  n2=MAF^2
  n1=2*MAF*(1-MAF)
  n0=(1-MAF)^2

  w2=1
  w1=n1/n2
  w0=n0/n2

  alpha=typeI/nTests

  if(verbose)
  {
    cat("mu2=", gm1, ", mu1=", gm2, ", mu0=", gm3, "\n")
    cat("n2=", n2, ", n1=", n1, ", n0=", n0, "\n")
    cat("alpha=", alpha, "\n")
  }

  k=3
  mydf1=k-1
  mydf2=myntotal-k
  q=qf(p=1-alpha, df1=mydf1, df2=mydf2)

  wVec=c(w2, w1, w0)
  wVec=wVec/sum(wVec, na.rm=TRUE)
  muVec=c(gm1, gm2, gm3)

  mu=sum(wVec*muVec, na.rm=TRUE)

  myncp = myntotal*sum(wVec*(muVec-mu)^2, na.rm=TRUE)
  myncp=myncp/(mystddev^2)

  power=1-pf(q=q, df1=mydf1, df2=mydf2,
             ncp=myncp)

  return(power)

}


# effsize=delta/sigma
powerEQTL.ANOVA2=function(effsize,
                          MAF,
                          typeI=0.05,
                          nTests=200000,
                          myntotal=200,
                          verbose=TRUE)
{

  myncp=2*MAF*(1-MAF)*myntotal*effsize^2

  alpha=typeI/nTests
  
  if(verbose)
  {
    cat("myncp=", myncp, "\n")
    cat("alpha=", alpha, "\n")
  }
  
  k=3
  mydf1=k-1
  mydf2=myntotal-k
  q=qf(p=1-alpha, df1=mydf1, df2=mydf2)
  
  power=1-pf(q=q, df1=mydf1, df2=mydf2,
             ncp=myncp)
  
  return(power)
  
}


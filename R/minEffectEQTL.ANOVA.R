# created on Dec. 11, 2016
# calculate minimum detectable effect size (delta/sigma) given sample size and power
#
# created on Dec. 8, 2016
#  (1) power calculation for eQTL analysis based on ANOVA or simple linear regression
#

diffPowerFunc.ANOVA2=function(effsize,
                              myntotal,
                              MAF=0.1,
                              typeI=0.05,
                              nTests=200000,
                              desiredPower=0.8)
{
  estPower=powerEQTL.ANOVA2(effsize=effsize,
                            MAF=MAF,
                           typeI=typeI,
                           nTests=nTests,
                           myntotal=myntotal,
                           verbose=FALSE)
  diff=(estPower-desiredPower)
  return(diff)

}


# MAF - minor allele frequency
# typeI - type I error rate
# nTests - number of tests
# myntotal - total number of subjects
# mystddev - standard deviation of gene expression levels
#            (assume each group of subjects has the same mystddev)
# deltaVec - mean difference of gene expression levels between groups
#            deltaVec[1]=mu2-mu1 and deltaVec[2]=mu3-mu2
# verbose - flag indicating if we should output intermediate results

minEffectEQTL.ANOVA=function(MAF, 
                   typeI=0.05,
                   nTests=200000,
                   myntotal=200,
                   mypower = 0.8,
                   verbose=TRUE)
{
  res.uniroot=uniroot(f=diffPowerFunc.ANOVA2,
                  interval=c(0.0001, 1000),
                  myntotal=myntotal,
                  typeI=typeI,
                  MAF = MAF,
                  nTests=nTests,  
                  desiredPower=mypower
                  )
  
  if(verbose)
  {
    cat("Results of uniroot>>>\n")
    print(res.uniroot)
  }
  effsize=res.uniroot$root
  
  return(effsize)
}


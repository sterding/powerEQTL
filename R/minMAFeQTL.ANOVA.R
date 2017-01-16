# created on Jan. 15, 2017
# calculate MAF based on simple linear regression
#

diffPowerFunc.ANOVA.MAF=function(MAF,
                              effsize,
                              myntotal,
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

minMAFeQTL.ANOVA=function(effsize, 
                   typeI=0.05,
                   nTests=200000,
                   myntotal=200,
                   mypower = 0.8,
                   verbose=TRUE)
{
  res.uniroot=uniroot(f=diffPowerFunc.ANOVA.MAF,
                  interval=c(0.000001, 0.5),
                  effsize = effsize,
                  myntotal=myntotal,
                  typeI=typeI,
                  nTests=nTests,  
                  desiredPower=mypower
                  )
  
  if(verbose)
  {
    cat("Results of uniroot>>>\n")
    print(res.uniroot)
  }
  MAF=res.uniroot$root
  
  return(MAF)
}


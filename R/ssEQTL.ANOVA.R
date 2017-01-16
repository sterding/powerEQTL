# created on Dec. 11, 2016
# (1) allow mu2-mu1 not equal to mu3-mu2,
#   where mu1 is the mean expression level for mutation homozygote,
# mu2 is the mean expression level for heterozygote,
#  and mu3 is the mean expression level for wildtype homozygote
#
# created on Dec. 8, 2016
#  (1) sample size calculation for eQTL analysis based on ANOVA or simple linear regression
#

# squared difference between estimated power and desired power
diffPowerFunc.ANOVA=function(myntotal,
                       MAF=0.1,
                       typeI=0.05,
                       nTests=200000,
                       mystddev=0.13,
                       deltaVec=c(0.13,0.13),
                       desiredPower=0.8)
{
  estPower=powerEQTL.ANOVA(MAF=MAF,
                              typeI=typeI,
                              nTests=nTests,
                              myntotal=myntotal,
                              mystddev=mystddev,
                              deltaVec=deltaVec, verbose=FALSE)
  diff=(estPower-desiredPower)^2
  return(diff)

}

# MAF - minor allele frequency
# typeI - type I error rate
# nTests - number of tests
# mypower - desired power
# mystddev - standard deviation of gene expression levels
#            (assume each group of subjects has the same mystddev)
# deltaVec - mean difference of gene expression levels between groups
#            deltaVec[1]=mu2-mu1 and deltaVec[2]=mu3-mu2
ssEQTL.ANOVA=function(MAF,
                   typeI=0.05,
                   nTests=200000,
                   mypower=0.8,
                   mystddev=0.13,
                   deltaVec=c(0.13,0.13))
{
  res.optim=optim(par=200, fn=diffPowerFunc.ANOVA,typeI=typeI,
                  nTests=nTests, mystddev=mystddev,
                  deltaVec=deltaVec, desiredPower=mypower,
        method = c("L-BFGS-B"),
        lower = 1, upper = 1000000)
  myntotal=res.optim$par

  return(myntotal)
}



diffPowerFunc.ANOVA3=function(
                              myntotal,
                              MAF=0.1,
                              typeI=0.05,
                              nTests=200000,
                              effsize=1,
                              desiredPower=0.8)
{
  estPower=powerEQTL.ANOVA2(effsize=effsize,
                            MAF=MAF,
                            typeI=typeI,
                            nTests=nTests,
                            myntotal=myntotal,
                            verbose=FALSE)
  diff=(estPower-desiredPower)^2
  return(diff)
  
}


ssEQTL.ANOVA2=function(
  effsize,
  MAF,
  typeI=0.05,
  nTests=200000,
  mypower=0.8
)
{
  res.optim=optim(par=200, fn=diffPowerFunc.ANOVA3,
                  MAF=MAF,
                  typeI=typeI,
                  nTests=nTests, 
                  effsize=effsize,
                  desiredPower=mypower,
                  method = c("L-BFGS-B"),
                  lower = 1, upper = 1000000)
  myntotal=res.optim$par
  
  return(myntotal)
}



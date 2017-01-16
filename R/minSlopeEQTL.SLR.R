# created on Dec. 8, 2016
#  (1) minimum slope calculation for eQTL analysis based on simple linear regression
#

# MAF - minor allele frequency
# typeI - type I error rate
# nTests - number of tests
# myntotal - total number of subjects
# mypower - desired power
# mystddev - standard deviation of gene expression levels
#            (assume each group of subjects has the same mystddev)
# verbose - flag indicating if intermedaite results should be output
minSlopeEQTL.SLR=function(MAF,
                   typeI=0.05,
                   nTests=200000,
                   myntotal=200,
                   mypower=0.8,
                   mystddev=0.13,
                   verbose=TRUE)
{

  sigma.x=sqrt(2*MAF*(1-MAF))
  alpha = typeI/nTests

  aa=minEffect.SLR(n=myntotal,
                   power=mypower,
                   sigma.x=sigma.x,
                   sigma.y=mystddev,
                   alpha = alpha,
                   verbose = verbose)
  return(aa$lambda.a)

}


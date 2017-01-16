# created on Dec. 8, 2016
#  (1) power calculation for eQTL analysis based on simple linear regression
#

# MAF - minor allele frequency
# typeI - type I error rate
# nTests - number of tests
# myntotal - total number of subjects
# mypower - desired power
# mystddev - standard deviation of gene expression levels
#            (assume each group of subjects has the same mystddev)

powerEQTL.SLR=function(MAF,
                          typeI=0.05,
                          nTests=200000,
                          slope=0.13,
                          myntotal=200,
                          mystddev=0.13,
                          verbose=TRUE)
{

  sigma.x=sqrt(2*MAF*(1-MAF))
  alpha = typeI/nTests

  aa=power.SLR(n=myntotal,
                   lambda.a=slope,
                   sigma.x=sigma.x,
                   sigma.y=mystddev,
                   alpha = alpha,
                   verbose = verbose)
  return(aa$power)

}


# created on Dec. 8, 2016
#  (1) sample size calculation for eQTL analysis based on simple linear regression
#

# MAF - minor allele frequency
# typeI - type I error rate
# nTests - number of tests
# myntotal - total number of subjects
# mypower - desired power
# mystddev - standard deviation of gene expression levels
#            (assume each group of subjects has the same mystddev)
# n.lower - lower bound for sample size
# n.upper - upper bound for sample size
# verbose - flag indicating if intermedaite results should be output
ssEQTL.SLR=function(MAF,
                       typeI=0.05,
                       nTests=200000,
                       slope=0.13,
                       mypower=0.8,
                       mystddev=0.13,
                       n.lower = 2.01,
                       n.upper = 1e+30,
                       verbose=TRUE)
{

  sigma.x=sqrt(2*MAF*(1-MAF))
  alpha = typeI/nTests

  aa=ss.SLR(power=mypower,
               lambda.a=slope,
               sigma.x=sigma.x,
               sigma.y=mystddev,
               n.lower = n.lower,
               n.upper = n.upper,
               alpha = alpha,
               verbose = verbose)
  return(aa$n)

}


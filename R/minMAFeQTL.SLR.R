# created on Jan. 15, 2017
#  (1) minimum MAF calculation for eQTL analysis based on simple linear regression
#



diffPowerFunc.SLR=function(MAF,
                           slope,   
                           myntotal,
                           mystddev=0.13,
                              typeI=0.05,
                              nTests=200000,
                              desiredPower=0.8)
{
  estPower=powerEQTL.SLR(MAF=MAF,
                            typeI=typeI,
                            nTests=nTests,
                            slope = slope,
                            myntotal=myntotal,
                            mystddev=mystddev,
                            verbose=FALSE)
  diff=(estPower-desiredPower)
  return(diff)
  
}




# slope - slope of the simple linear regression
# typeI - type I error rate
# nTests - number of tests
# myntotal - total number of subjects
# mypower - desired power
# mystddev - standard deviation of gene expression levels
#            (assume each group of subjects has the same mystddev)
# verbose - flag indicating if intermedaite results should be output
minMAFeQTL.SLR=function(slope,
                   typeI=0.05,
                   nTests=200000,
                   myntotal=200,
                   mypower=0.8,
                   mystddev=0.13,
                   verbose=TRUE)
{

  res.uniroot=uniroot(f=diffPowerFunc.SLR,
    interval=c(0.000001, 0.5),
    slope=slope,   
    myntotal=myntotal,
    mystddev=mystddev,
    typeI=typeI,
    nTests=nTests,
    desiredPower=mypower)
                        
  if(verbose)
  {
    cat("Results of uniroot>>>\n")
    print(res.uniroot)
  }
  MAF=res.uniroot$root
  
  return(MAF)
}


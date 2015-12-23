# This function defines and runs the model
# The data is not standardized
genMCMC = function( data, yName="y" , 
                    numSavedSteps=50000 , saveName=NULL ) { 
  require(rjags)
  #-----------------------------------------------------------------------------
  # THE DATA.
  y = data[,yName]
  #x = data[,xName]
  # Specify the data in a list, for later shipment to JAGS:
  dataList = list(
   # x = x ,
    y = y 
  )
  #-----------------------------------------------------------------------------
  # THE MODEL.
  # Specify the model for unstandardized data:
  modelString = "
  data {
    Ntotal <- length(y)
   # xm <- mean(x)
    #ym <- mean(y)
  }
  model {
  for ( i in 1:Ntotal ) {
  y[i] ~ dbern( theta )
  }
  # Priors vague 
  theta ~ dbeta(1,1)  
  }
  " # close quote for modelString
  # Write out modelString to a text file
  writeLines( modelString , con="TEMPmodel.txt" )
  #-----------------------------------------------------------------------------
  # INTIALIZE THE CHAINS.
  # Let JAGS do it...
  #-----------------------------------------------------------------------------
  # RUN THE CHAINS
  parameters = c( "theta" )
  adaptSteps = 500  # Number of steps to "tune" the samplers
  burnInSteps = 1000
  nChains = 4 
  thinSteps = 1
  nIter = 5000
  #nIter = ceiling( ( numSavedSteps * thinSteps ) / nChains )
  # Create, initialize, and adapt the model:
  jagsModel = jags.model( "TEMPmodel.txt" , data=dataList , #inits=initsList , 
                          n.chains=nChains , n.adapt=adaptSteps )
  # Burn-in:
  cat( "Burning in the MCMC chain...\n" )
  update( jagsModel , n.iter=burnInSteps )
  # The saved MCMC chain:
  cat( "Sampling final MCMC chain...\n" )
  codaSamples = coda.samples( jagsModel , variable.names=parameters , 
                              n.iter=nIter , thin=thinSteps )
  # resulting codaSamples object has these indices: 
  #   codaSamples[[ chainIdx ]][ stepIdx , paramIdx ]
  if ( !is.null(saveName) ) {
    save( codaSamples , file=paste(saveName,"Mcmc.Rdata",sep="") )
  }
  return( codaSamples )
} # end function
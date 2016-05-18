# This function defines and runs the model
# The data is not standardized
genMCMC = function( data ,  yName="y" , sName="s"  , 
                    numSavedSteps=50000 , saveName=NULL ) { 
  require(rjags)
  #-----------------------------------------------------------------------------
  # THE DATA.
  y = data[,yName]
  s = as.numeric(factor(data[,sName]))
 # x = data[,xName]
  # Specify the data in a list, for later shipment to JAGS:
  dataList = list(
   # x = x ,
    y = y, 
    s = s,
    Ntotal = length(y),
    Nsubj = length(unique(s))
  )
  #-----------------------------------------------------------------------------
  # THE MODEL.
  # Specify the model for unstandardized data:
  modelString = "
  model {
  for ( i in 1:Ntotal ) {
  y[i] ~ dt( beta[s[i]], 1/sigma^2, 4 )
  }
  for ( s in 1:Nsubj ) {
      beta[s] ~ dnorm( betamu , 1/(betasigma)^2 )  
    }
  # Priors vague 
  betamu ~ dnorm( 0 , 1/(10)^2 )
  betasigma ~ dunif( 1.0E-3 , 1.0E+2 )
  sigma ~ dunif( 1.0E-3 , 1.0E+3 )
  }
  " # close quote for modelString
  # Write out modelString to a text file
  writeLines( modelString , con="TEMPmodel.txt" )
  #-----------------------------------------------------------------------------
  # INTIALIZE THE CHAINS.
  # Let JAGS do it...
  #-----------------------------------------------------------------------------
  # RUN THE CHAINS
  #parameters = c( "beta0" ,  "beta1" ,  "sigma"  )
  parameters = c( "beta" ,   "sigma"  )
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
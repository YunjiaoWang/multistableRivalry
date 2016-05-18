# This function defines and runs the model
# The data is not standardized
genMCMC_trans = function( data,  xName="x", yName="y", sName = "s",
                    numSavedSteps=50000 , saveName=NULL ) { 
  require(rjags)
  #-----------------------------------------------------------------------------
  # THE DATA.
  y = data[,yName]
  x = data[,xName]
  s = data[,sName]
  
  # Specify the data in a list, for later shipment to JAGS:
  dataList = list(
    x = x ,
    y = y ,
    s = s ,
    Ntotal = length(y),
    Nsubj = length(unique(s))
  )
  #-----------------------------------------------------------------------------
  # THE MODEL.
  # Specify the model for unstandardized data:
  modelString = "
  model {
  for ( i in 1:Ntotal ) {
  y[i] ~ dbern( theta0[s[i]] + theta1[s[i]]*x[i] )
  }
  for ( s in 1:Nsubj ) {
      # Use the mode (omega) and concentration (kappa) parmaters for the beta
      theta0[s] ~ dbeta(omega*(kappa-2)+1,(1-omega)*(kappa-2)+1)  
      theta1[s] ~ dnorm( theta1mu , 1/(theta1sigma)^2 )
    }
  # Priors vague 
  theta1mu ~ dnorm( 0 , 1/(10)^2 )
  theta1sigma ~ dunif( 1.0E-3 , 1.0E+2 )
  omega ~ dbeta(1,1)
  kappa <- kappaMinusTwo + 2
  kappaMinusTwo ~ dgamma(0.01,0.01)
  }
  " # close quote for modelString
  # Write out modelString to a text file
  writeLines( modelString , con="TEMPmodel.txt" )
  #-----------------------------------------------------------------------------
  # INTIALIZE THE CHAINS.
  # Let JAGS do it...
  #-----------------------------------------------------------------------------
  # RUN THE CHAINS
  parameters = c( "theta0","theta1" )
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
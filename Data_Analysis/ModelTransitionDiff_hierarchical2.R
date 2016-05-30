# This function defines and runs the model
# The data is not standardized
genMCMC_trans = function( data , xName="x" , yName="y" , sName="s" ,
                          numSavedSteps=50000 , thinSteps = 1 , 
                          saveName=NULL , runjagsMethod=runjagsMethodDefault ,
                          nChains=4) 
{
  require(rjags)
  #-----------------------------------------------------------------------------
  # THE DATA.
  y = data[,yName]
  x = data[,xName]
  s = data[,sName]
  Ntotal = length(y)
  Nsubj = length(unique(s))
  
  # Specify the data in a list, for later shipment to JAGS:
  dataList = list(
    x = x ,
    y = y ,
    s = s ,
    Ntotal = Ntotal,
    Nsubj = Nsubj
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
  parameters = c( "theta0","theta1","theta1mu","theta1sigma","omega","kappa" )
  
  adaptSteps = 5000  # Number of steps to "tune" the samplers
  burnInSteps = 5000
  runJagsOut <- run.jags( method=runjagsMethod ,
                          model="TEMPmodel.txt" , 
                          monitor=parameters , 
                          data=dataList ,  
                          #inits=initsList , 
                          n.chains=nChains ,
                          adapt=adaptSteps ,
                          burnin=burnInSteps , 
                          sample=ceiling(numSavedSteps/nChains) ,
                          thin=thinSteps ,
                          summarise=FALSE ,
                          plots=FALSE )
  codaSamples = as.mcmc.list( runJagsOut )
  
  return( codaSamples )
} # end function

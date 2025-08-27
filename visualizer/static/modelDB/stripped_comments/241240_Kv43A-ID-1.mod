NEURON {
  SUFFIX 	Kv43A 
  USEION 	k READ ek WRITE ik
  RANGE		gbar, g, i

  GLOBAL 	inf_n, tau_n, inf_h, tau_h
}

PARAMETER {
  ek 			(mV)
  gbar		= 30 	(pS/um2)

  
  gates_n	= 4		
  vhalf_n	= -32	(mV)
  slope_n	= -6.0	(mV)	
  needAdj = 1
  vhalfA_n 	= 0 (mV)		
  slopeA_n	= 0 (mV)
  v5_adj	= 0 (mV)		
  slp_adj	= 0 (mV)

  
  tauA_n	= 2.5	(ms)
  tauDv_n	= 0	(mV)	
  tauG_n	= 0.4		
  tauF_n	= 0		
  tau0_n	= 0.2	(ms)	

 
  vhalf_h	= -64  (mV)
  slope_h	= 6.5	(mV)	

 
  tauA_h	= 20	(ms)
  tauDv_h	= 0	(mV)	
  tauG_h	= 0.5		
  tauF_h	= 0		
  tau0_h	= 4.0	(ms)	
}

STATE {
  n
  h
}

ASSIGNED {
  v		(mV)
  celsius	(degC)
  ik 		(mA/cm2)
  g		(pS/um2)
  i		(mA/cm2)
  inf_n
  tau_n		(ms)
  inf_h
  tau_h		(ms)
}

BREAKPOINT {
  SOLVE states METHOD cnexp
  g = 0
  if( n >=0 ){	
    g 	= gbar * n^gates_n * h
  }
  i	= g * ( v - ek ) * (1e-4)
  ik	= i
}

INITIAL {
  needAdj = 1
  rates( v )
  n = inf_n
  h = inf_h
}

UNITS {
  (mA)	= (milliamp)
  (mV)	= (millivolt)
  (pS)	= (picosiemens)
  (um)	= (micrometer)
}

DERIVATIVE states {     
  rates( v )
  n' = ( inf_n - n )/ tau_n
  h' = ( inf_h - h )/ tau_h
}

PROCEDURE rates( v (mV)){
  if( needAdj > 0 ){
    needAdj = 0
    ngate_adjust( gates_n, vhalf_n, slope_n )
    vhalfA_n = v5_adj
    slopeA_n = slp_adj
  }
  inf_n = Boltzmann( v, vhalfA_n, slopeA_n )
  tau_n = BorgMod_tau( v, vhalfA_n, slopeA_n, tau0_n, tauA_n, tauG_n, tauF_n, tauDv_n )
  
  inf_h = Boltzmann( v, vhalf_h, slope_h )
  tau_h = BorgMod_tau( v, vhalf_h, slope_h, tau0_h, tauA_h, tauG_h, tauF_h, tauDv_h )
}

FUNCTION Boltzmann( v (mV), v5 (mV), s (mV) ){
  Boltzmann = 1 / (1 + exp( (v - v5) / s ))
}

FUNCTION BorgMod_tau( v (mV), v5 (mV), s (mV), tau0 (ms), tauA (ms), tauG, tauF, tauDv (mV) ) (ms) {
  LOCAL kc, kr, Dv, wr, kf



  Dv = (v - ( v5 + tauDv ))
  kf =  10^tauF


  BorgMod_tau = tau0 + tauA * 4 * sqrt( tauG * (1-tauG))
 	        / ( exp( - Dv *(1-tauG)*kf/s ) + exp( Dv *tauG*kf/s ))
}


FUNCTION Boltz_m1( x, v5 (mV), s (mV) ) (mV) {
  Boltz_m1 = s * log( 1/x - 1 ) + v5
}




PROCEDURE ngate_adjust( ng, vh (mV), slp (mV) ) {
  LOCAL x1, x2, v1, v2
  x1 = 0.3
  x2 = 0.7
  v1 = Boltz_m1( x1, vh, slp )
  v2 = Boltz_m1( x2, vh, slp )
  slp_adj = (v2 - v1)/( log( (1/x2)^(1/ng) - 1 ) - log( (1/x1)^(1/ng) - 1 ) )
  v5_adj = v1 - slp_adj * log( 1 / x1^(1/ng) - 1 )
}
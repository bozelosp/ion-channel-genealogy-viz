NEURON {
  SUFFIX 	H_HCNl
  NONSPECIFIC_CURRENT i
  RANGE		gbar, g, eh

  GLOBAL 	inf_n, tau_n, inf_ns, tau_ns
}

PARAMETER {
  eh		= -32	(mV)
  gbar		= 30 	(pS/um2)

  
  gates_n	= 2
  vhalf_n	= -98	(mV)
  slope_n	= 9	(mV)
  needAdj = 1
  vhalfA_n 	= 0 (mV)		
  slopeA_n	= 0 (mV)
  v5_adj	= 0 (mV)		
  slp_adj	= 0 (mV)



  
  gates_ns	= 1
  vhalf_ns	= -98	(mV)
  slope_ns	= 9	(mV)







  
  tauA_n	= 107	(ms)
  tauDv_n	= 0	(mV)	
  tauG_n	= 0.37		
  tauF_n	= 0		
  tau0_n	= 14.7 	(ms)	

  
  tauA_ns	= 367	(ms)
  tauDv_ns	= 0	(mV)	
  tauG_ns	= 0.3975	
  tauF_ns	= 0		
  tau0_ns	= 161.3	 (ms)	

}

STATE {
  n
  ns
}

ASSIGNED {
  v		(mV)
  celsius	(degC)
  ik 		(mA/cm2)
  g		(pS/um2)
  i		(mA/cm2)
  inf_n
  tau_n		(ms)
  inf_ns
  tau_ns	(ms)
}

BREAKPOINT {
  SOLVE states METHOD cnexp
  g = 0
  if( n >=0 ){	
    g 	= gbar * (0.434*n^gates_n+(1-0.434)*ns^gates_ns)
  }
  i	= g * ( v - eh ) * (1e-4)
}

INITIAL {
  needAdj = 1
  rates( v )
  n = inf_n
  ns = inf_ns
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
  ns' = ( inf_ns - ns )/ tau_ns
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
  inf_ns = Boltzmann( v, vhalf_ns, slope_ns )
  tau_ns = BorgMod_tau( v, vhalf_ns, slope_ns, tau0_ns, tauA_ns, tauG_ns, tauF_ns, tauDv_ns )
  
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
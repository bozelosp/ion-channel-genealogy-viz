NEURON {
  POINT_PROCESS GABA
  RANGE gbar, tau_r, tau_d, scale, spkcnt, countflag, i, t1, Erev, qfact
  NONSPECIFIC_CURRENT i
}

UNITS {
  (nA) = (nanoamp)
  (mV) = (millivolt)
  (umho) = (micromho)
}

PARAMETER {
	gbar = 0.0021  (umho)	
							
	tau_r = 0.5 	(ms)   	
	tau_d = 7.5  	(ms)   	
	Erev  = -60    (mV)   	
	saturate = 1.2 			
							
	qfact = 2				
}							


ASSIGNED {
	g (umho)
	v (mV)   		
	i (nA)   		
	t1  (ms)
	
	y1_add (/ms)    
	y1_loc (/ms)

	countflag	
	spkcnt		

	scale		
}			

STATE { 
	y1 (/ms) 
	y2    			
}

INITIAL {
  y1_add = 0
	scale = 0
	spkcnt = 0
	countflag = 0
	t1 = 0
	y1_loc = 0
}

BREAKPOINT {
  SOLVE betadyn METHOD cnexp
	g =  gbar * y2 
  i = scale * g * (v - Erev)
}

DERIVATIVE betadyn {
  
  y1' = -y1 / (tau_d/qfact)
  y2' = y1 - y2 / (tau_r/qfact)
}

NET_RECEIVE( weight, y1_loc (/ms)) {
  
  y1_loc = y1_loc*exp( -(t - t1) / (tau_d/qfact) )

  
  
  y1_add = (1 - y1_loc/saturate)

  
  y1_loc = y1_loc + y1_add

  
  y1 = y1 + y1_add

  
  t1 = t

	spkcnt = spkcnt + 1

	scale = weight
}
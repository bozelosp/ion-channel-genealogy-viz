NEURON {
	POINT_PROCESS AMPA
	RANGE gbar, tau_r, tau_d, scale, spkcnt, countflag, i, t1, ca_ratio, ical, itmp, qfact, g
	NONSPECIFIC_CURRENT i
 	USEION cal WRITE ical VALENCE 2

}

UNITS {
	(nA) = (nanoamp)
	(mV) = (millivolt)
	(umho) = (micromho)
}

PARAMETER {
	gbar = 8.5e-4   (umho) 	
							
							
	tau_r = 2.2 	(ms)   	
	tau_d = 11.5  	(ms)   	
	
	Erev = 0    	(mV)   	
	saturate = 1.2 			
							
	qfact = 2				
	ca_ratio = 0.005		
                            
							
    g_factor                

}							
	


ASSIGNED {
	g (umho)
	v (mV)   		
	itmp	(nA)	
	i (nA)   		
	ical (nA)		
	t1 (ms)
	
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
	g_factor = 1
}

BREAKPOINT {
  	SOLVE betadyn METHOD cnexp
	g = gbar * g_factor * y2
  	itmp = scale * g * (v - Erev)
  	i = (1-ca_ratio) * itmp
  	ical = ca_ratio * itmp
}

DERIVATIVE betadyn {
	
	y1' = -y1 / (tau_d/qfact)
	y2' = y1 - y2 / (tau_r/qfact)
}

NET_RECEIVE( weight, y1_loc (/ms) ) {
	
	y1_loc = y1_loc*exp( -(t - t1) / (tau_d/qfact) )

	
	
	y1_add = (1 - y1_loc/saturate)

	
	y1_loc = y1_loc + y1_add

	
	y1 = y1 + y1_add

	
	t1 = t

	spkcnt = spkcnt + 1

	scale = weight
}
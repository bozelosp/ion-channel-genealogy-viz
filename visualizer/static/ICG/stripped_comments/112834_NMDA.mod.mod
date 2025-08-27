NEURON {
	POINT_PROCESS NMDA
	RANGE gbar, ca_ratio, tau_r, tau_d, scale, spkcnt, countflag, mg, i, ical, t1, itmp, qfact
	NONSPECIFIC_CURRENT i
	USEION cal WRITE ical VALENCE 2
}

UNITS {
	(nA) = (nanoamp)
	(mV) = (millivolt)
	(umho) = (micromho)
	(mM) = (milli/liter)
}

PARAMETER {
	gbar = 12.2e-5 (umho)
						
	tau_r = 5.63 (ms)   
	tau_d = 320  (ms)   
	
	Erev  = 0    (mV)   
	saturation = 7.0 	
						
	qfact = 2			
	mg = 1      (mM)    
	ca_ratio = 0.1		
						
}

ASSIGNED {
	g (umho)
	v (mV)   		
	itmp	(nA)	
	i (nA)   		
	ical (nA)		
	t1	(ms)
	
	y1_add (/ms)    
	y1_loc (/ms)

	spkcnt		
	countflag	

	B				
	scale		
}			


STATE { 
	y1 (/ms) 
	y2    			
}

INITIAL {
	y1_add = 0
  
  	B = mgblock(v)
	scale = 1
	spkcnt = 0
	countflag = 0
	t1 = 0
	y1_loc = 0
}

BREAKPOINT {
	SOLVE betadyn METHOD cnexp
  	mgblock(v)
	g = gbar * y2
	itmp = scale * g * B * (v - Erev)	
	i = (1-ca_ratio) * itmp
	ical = ca_ratio * itmp
}

DERIVATIVE betadyn {
  
	y1' = -y1 / (tau_d/qfact)
	y2' = y1 - y2 / (tau_r/qfact)
}

NET_RECEIVE( weight, y1_loc (/ms)) {
	
	y1_loc = y1_loc*exp( -(t - t1) / (tau_d/qfact) )

	
	
	y1_add = (1 - y1_loc/saturation)

	
	y1_loc = y1_loc + y1_add

	
	y1 = y1 + y1_add

	
	t1 = t

	spkcnt = spkcnt + 1

	scale = weight
}


PROCEDURE mgblock( v(mV) ) {
	

	TABLE B DEPEND mg
		FROM -100 TO 100 WITH 201

	B = 1 / (1 + exp(0.062 (/mV) * -v) * (mg / 3.57 (mM)))
}
NEURON {
	
	SUFFIX nap
	USEION na READ ena WRITE ina
	RANGE  gbar, timestauh, timestaum, shifttaum, shifttauh
	GLOBAL minf, mtau
	THREADSAFE
}

PARAMETER {
	
	gbar = 5   	(pS/um2)
	
	
	
	
	timestauh=1
	timestaum=1
	shifttaum=1
	shifttauh=1
	
	
	ena	= 55	(mV)		
	
	celsius 	(degC)
	v 			(mV)
}


UNITS {
	
	(mA) = (milliamp)
	(mV) = (millivolt)
	(pS) = (picosiemens)
	(um) = (micron)
} 

ASSIGNED {
	
	ina 		(mA/cm2)
	inap		(mA/cm2)
	minf 		
	mtau (ms)	
}

STATE { m }

UNITSOFF

BREAKPOINT {
    
    SOLVE states METHOD cnexp
	mtau = 1
	minf = (1/(1+exp(-(v+62.3)/6.8))) 
	ina = (1e-4)*gbar*m * (v - ena)
} 

INITIAL {
	
	mtau = 1
	minf = (1/(1+exp(-(v+62.3)/6.8))) 
	m=minf  
}

DERIVATIVE states {   
  
	mtau = 1
	minf = (1/(1+exp(-(v+62.3)/6.8))) 
	m' = (minf-m)/mtau
}

UNITSON
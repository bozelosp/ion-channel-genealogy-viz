NEURON {
	SUFFIX naIn
	USEION na READ ena WRITE ina
	RANGE  gbar, thegna, htau
	GLOBAL minf, hinf 
}

PARAMETER {
	gbar = .0052085   	(mho/cm2)
	
	
	
	
	mtau = 1 (ms)
	
	
	eNa = 55 	(mV)		
	ena		(mV)            
	celsius (degC)
	v 		(mV)
}


UNITS {
	(mA) = (milliamp)
	(mV) = (millivolt)
	(pS) = (picosiemens)
	(um) = (micron)
} 

ASSIGNED {
	ina 		(mA/cm2)
	thegna		(mho/cm2)
	minf 		hinf 		
	htau (ms)	

}
 

STATE { m h}



BREAKPOINT {
        SOLVE states METHOD cnexp
	trates(v)	
	
	thegna =gbar*m*h       
	ina = thegna * (v - eNa)
	} 

INITIAL {
	trates(v)
	m=minf  
	h=hinf
}

DERIVATIVE states {   
        
	m' = (minf-m)/mtau
        h' = (hinf-h)/htau

}



FUNCTION alphah(vm (mV)) (/ms/mV)  { LOCAL a, b, k
  UNITSOFF
  a = -0.00288
  b = -0.049
  k = 4.63	
  alphah = 1000*(a*vm+b)/(1-exp((vm+b/a)/k)) 
  UNITSON
}



FUNCTION betah(vm (mV)) (/ms/mV) { LOCAL a, b, k
  UNITSOFF
  a = 0.00694
  b = 0.447
  k = -2.63
  betah = 1000*(a*vm+b)/(1-exp((vm+b/a)/k)) 
  UNITSON
}

PROCEDURE trates(vm (mV)) {LOCAL alpha, beta
	  
	UNITSOFF
        minf = (1/(1+exp(-(v+52.6)/4.6))) 
	alpha = alphah(vm)
	beta = betah(vm)	
	hinf = alpha/(alpha+beta) 
	htau = 1/(alpha+beta)	
	UNITSON
}


UNITSON
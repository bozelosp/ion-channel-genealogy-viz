NEURON {
	SUFFIX nap
	USEION na READ ena WRITE ina
	RANGE  gbar, timestauh, timestaum, shifttaum, shifttauh, thegna
	GLOBAL minf, mtau 

}

PARAMETER {
	gbar = .0052085   	(mho/cm2)
	
	
	
	
	timestauh=1
	timestaum=1
	shifttaum=1
	shifttauh=1
	
	
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
	minf 		
	mtau (ms)	
}
 

STATE { m }


UNITSOFF

BREAKPOINT {
    SOLVE states METHOD cnexp
	mtau = 1
	minf = (1/(1+exp(-(v+52.3)/6.8))) 
	
	thegna =gbar*m       

	ina = thegna * (v - eNa)
	} 

INITIAL {

	mtau = 1
	minf = (1/(1+exp(-(v+52.3)/6.8))) 
	m=minf  

}

DERIVATIVE states {   
  
	mtau = 1
	minf = (1/(1+exp(-(v+52.3)/6.8))) 
	m' = (minf-m)/mtau


}



UNITSON
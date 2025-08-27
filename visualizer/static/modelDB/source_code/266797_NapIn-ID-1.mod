TITLE nap
:persisten sodium current with inactivation

NEURON {
	SUFFIX napIn
	USEION na READ ena WRITE ina
	RANGE  gbar, thegna, htau
	GLOBAL minf, mtau, hinf, mintaua, mintaub
}

PARAMETER {
	gbar = .0052085   	(mho/cm2)
    htau = 4000	(ms)
	
	eNa = 55 	(mV)		
	ena		(mV)            
	celsius (degC)
	v 		(mV)
	sh0=6
	sh1=0
	km=2
	kl=10
	mintaua=1000
	mintaub=100
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
	mtau (ms)		
}
 

STATE { m h}

: hier eigener Befehl
UNITSOFF

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
    trates(v)      
	if (m<minf) {mtau=mintaua} else {mtau=mintaub}
	m' = (minf-m)/mtau
    h' = (hinf-h)/htau

}

PROCEDURE trates(vm) {  
	minf = 1/(1+exp(-(v+sh0+52.3)/km))
    hinf = 1/(1+exp((v+sh1+48)/kl))
	
}



UNITSON

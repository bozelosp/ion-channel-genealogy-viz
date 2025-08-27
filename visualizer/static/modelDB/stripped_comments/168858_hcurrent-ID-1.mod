UNITS { 
	(mV) = (millivolt) 
	(mA) = (milliamp) 
} 

NEURON { 
	SUFFIX ar
	NONSPECIFIC_CURRENT i
	RANGE gbar              
        RANGE Vh, tauMin
        RANGE a, b
}

PARAMETER { 
	gbar = 0.0 	(mho/cm2)
	erev = -55	(mV)  
        Vh = -75	(mV)
        a = -0.4090909	(/mV)
        b = 0.001	(1)
        tauMin = 5.0	(ms)
} 

ASSIGNED { 
	i 		(mA/cm2) 
	minf 		(1)
	mtau 		(ms) 
	v		(mV)
} 
STATE {
	m
}

INITIAL {
	rates(v)
	m = minf
}

BREAKPOINT {
	SOLVE states METHOD cnexp
	i = gbar * m * ( v - erev )
}

DERIVATIVE states {
	rates(v)
	m' = (minf-m)/mtau
}

UNITSOFF 

PROCEDURE rates(V (mV)) {
	minf = 1 / ( 1 + exp( -2 * a * ( V - Vh )) )
	mtau = 1 / b / (exp( a*(V-Vh)) + exp (-a*(V-Vh)) )

        if( mtau < tauMin ) { mtau = tauMin }
}

UNITSON
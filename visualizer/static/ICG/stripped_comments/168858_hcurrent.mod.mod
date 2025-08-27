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
	GLOBAL eh
}

PARAMETER { 
	gbar = 1.0 	(mho/cm2)
	
        Vh = -75	(mV)
        a = -0.4090909	(/mV)
        b = 0.001	(1)
        tauMin = 5.0	(ms)
} 

ASSIGNED {
	eh (mV) 
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
	i = gbar * m * ( v - eh )
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
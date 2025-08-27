INDEPENDENT { t FROM 0 TO 1 WITH 1 (ms) }

UNITS { 
	(mV) = (millivolt) 
	(mA) = (milliamp) 
	(mM) = (milli/liter)
}
 
NEURON { 
	SUFFIX kc_fast
	USEION k READ ek WRITE ik
	USEION ca READ cai
	RANGE  gbar, ik
}

PARAMETER { 
	gbar = 0.0 	(mho/cm2)
	v (mV) ek 		(mV)  
	cai		(mM)
} 

ASSIGNED { 
	ik 		(mA/cm2) 
	alpha (/ms) beta	(/ms)
}
 
STATE {
	m
}

BREAKPOINT { 
	SOLVE states METHOD cnexp
	if( 0.004(/mM) * cai < 1 ) {
		ik = gbar * m * 0.004(/mM) * cai * ( v - ek ) 
	}else{
		ik = gbar * m * ( v - ek ) 
	}
}
 
INITIAL { 
	settables(v) 
	m = alpha / ( alpha + beta )
	m = 0
}
 
DERIVATIVE states { 
	settables(v) 
	m' = alpha * ( 1 - m ) - beta * m 
}

UNITSOFF 

PROCEDURE settables(v(mV)) { 
	TABLE alpha, beta FROM -120 TO 40 WITH 641

	if( v <= -10.0 ) {
		alpha = 2 / 37.95 * ( exp( ( v + 50 ) / 11 - ( v + 53.5 ) / 27 ) )

		
		beta  = 2 * exp( ( - v - 53.5 ) / 27 ) - alpha
	}else{
		alpha = 2 * exp( ( - v - 53.5 ) / 27 )
		beta  = 0
	}
	
	alpha = alpha * 2
	beta  = beta  * 2
}

UNITSON
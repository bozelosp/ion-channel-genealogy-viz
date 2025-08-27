INDEPENDENT { t FROM 0 TO 1 WITH 1 (ms) }

NEURON { 
	SUFFIX iCat2
	USEION ca READ eca WRITE ica
	RANGE  gbar, ica, eca
}


UNITS { 
	(mV) = (millivolt) 
	(mA) = (milliamp) 
}
 

PARAMETER { 
	v  		(mV)  
	gbar = 1.0 	(mho/cm2)
	cai	= 8e-5 (mM)		
	cao	= 2    (mM)
}
 

STATE {
	m
}

ASSIGNED { 
	ica 		(mA/cm2) 
	alpha beta	(/ms)
	eca	(mV)
	i	(mA/cm2)
}
 


BREAKPOINT { 
	SOLVE states METHOD cnexp
	


	

	ica = gbar * m * m * ( v - eca) 
        i = ica		
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

PROCEDURE settables(v) { LOCAL tmp
	TABLE alpha, beta FROM -120 TO 40 WITH 641

	alpha = 1.11 / ( 1 + exp( - 0.058 * ( v +10 ) ) )
	tmp = v +23.9
	if ( fabs( tmp ) < 1e-6 ) {
		beta  = 0.1 * exp( - tmp / 5 ) 
	}else{
		beta  = 0.02 * tmp / ( exp( tmp / 5 ) - 1 )
	}
}

UNITSON
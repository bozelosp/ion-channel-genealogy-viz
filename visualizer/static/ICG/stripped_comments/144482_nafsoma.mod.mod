INDEPENDENT { t FROM 0 TO 1 WITH 1 (ms) }

UNITS { 
	(mV) = (millivolt) 
	(mA) = (milliamp) 
} 
NEURON { 
	SUFFIX nafsoma
	USEION na READ ena WRITE ina
	RANGE gbar, ina
	RANGE v11, v22, v33, v44,fastNashift
}
PARAMETER { 
	fastNashift = 0
	gbar = 1.0 	   (mho/cm2)
	v ena 		   (mV) 

	v11 = 62.9 (mV)
	v22 = 10.7 (mV)
	v33 = 37   (mV)
	v44	= 10   (mV)
} 
ASSIGNED { 
	ina 		   (mA/cm2) 
	minf hinf 	   (1)
	mtau htau 	   (ms) 
} 
STATE {
	m h
}
BREAKPOINT { 
	SOLVE states METHOD cnexp
	ina = gbar * m * m * m * h * ( v - ena ) 
} 
INITIAL { 
	settables( v - fastNashift ) 
	m = minf
	m = 0
	h  = hinf
} 
DERIVATIVE states { 
	
	
	minf  = 1 / ( 1 + exp( ( - ( v + fastNashift ) - 30 ) / 10 ) ) 
	if( ( v + fastNashift ) < -28.0 ) {
		mtau = 0.025 + 0.14 * exp( ( ( v + fastNashift ) + 28 ) / 10 )
	} else{
		mtau = 0.02 + 0.145 * exp( ( - ( v + fastNashift ) - 30 ) / 10 ) 
	}
	hinf  = 1 / ( 1 + exp( ( ( v + fastNashift  ) + v11 ) / v22 ) ) 
	htau = 0.15 + 1.15 / ( 1 + exp( ( ( v + fastNashift )+ v33 ) / v44 ) ) 
	
	
	m' = ( minf - m ) / mtau 
	h' = ( hinf - h ) / htau
}

UNITSOFF 

PROCEDURE settables(v1(mV)) {

	TABLE minf, hinf, mtau, htau  FROM -120 TO 40 WITH 641

	minf  = 1 / ( 1 + exp( ( - ( v1 + fastNashift ) - 30 ) / 10 ) ) 
	if( ( v1 + fastNashift ) < -28.0 ) {
		mtau = 0.025 + 0.14 * exp( ( ( v1 + fastNashift ) + 28 ) / 10 )
	} else{
		mtau = 0.02 + 0.145 * exp( ( - ( v1 + fastNashift ) - 30 ) / 10 ) 
	}

	









	
hinf  = 1 / ( 1 + exp( ( ( v1 + fastNashift  ) + v11 ) / v22 ) ) 
htau = 0.15 + 1.15 / ( 1 + exp( ( ( v1 + fastNashift )+ v33 ) / v44 ) ) 
	
	
}

UNITSON
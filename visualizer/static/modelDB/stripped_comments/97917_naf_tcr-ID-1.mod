INDEPENDENT { t FROM 0 TO 1 WITH 1 (ms) }

UNITS { 
	(mV) = (millivolt) 
	(mA) = (milliamp) 
} 
NEURON { 
	SUFFIX naf_tcr
	USEION na READ ena WRITE ina
	RANGE gbar, ina,m, h, df, shift_mnaf, minf, mtau
	RANGE shift_hnaf, shift_mnaf_init, shift_mnaf_run, hinf, htau
	RANGE shift_hnaf_run
}
PARAMETER { 
	shift_mnaf_init =-3 (mV) 
	shift_mnaf_run = -2.5 (mV) 
	shift_hnaf = -7.0 (mV)
	gbar = 0.0 	   (mho/cm2)
	v (mV) ena 		   (mV)  
} 
ASSIGNED { 
	shift_mnaf (mV)
	ina 		   (mA/cm2) 
	minf (1)
	hinf (1)
	mtau (ms)
	htau (ms)
	df (mV)
} 
STATE {
	m h
}
BREAKPOINT { 
	SOLVE states METHOD cnexp
	ina = gbar * m * m * m * h * ( v - ena ) 
	df = v - ena
}
INITIAL { 
	shift_mnaf = shift_mnaf_init + shift_mnaf_run
	settables( v )
	m = minf
	m = 0
	h  = hinf
} 
DERIVATIVE states { 
	settables( v ) 
	m' = ( minf - m ) / mtau 
	h' = ( hinf - h ) / htau
}

UNITSOFF 

PROCEDURE settables(v1(mV)) {

	TABLE minf, hinf, mtau, htau  FROM -120 TO 40 WITH 641

	minf  = 1 / ( 1 + exp( ( - ( v1 + shift_mnaf ) - 38 ) / 10 ) )
	if( ( v1 + shift_mnaf ) < -30.0 ) {
		mtau = 0.025 + 0.14 * exp( ( ( v1 + shift_mnaf ) + 30 ) / 10 )
	} else{
		mtau = 0.02 + (0.145) * exp( ( - ( v1 + shift_mnaf ) - 30 ) / 10 ) 
	}

	

	hinf  = 1.0 / ( 1.0 + exp( ( v1 + shift_hnaf + 62.9 ) / 10.7 ) )
	htau = 0.15 + 1.15 / ( 1.0 + exp( ( v1 + 37.0 ) / 15.0 ) )
}

UNITSON
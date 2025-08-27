INDEPENDENT { t FROM 0 TO 1 WITH 1 (ms) }

UNITS { 
	(mV) = (millivolt) 
	(mA) = (milliamp) 
}
 
NEURON { 
	SUFFIX cat_a
	
	USEION ca READ eca WRITE ica
        RANGE gbar, m, h, alphah, betah 	
}

PARAMETER { 
	gbar = 1.0 	(mho/cm2)
	v 		(mV)  
}
 
ASSIGNED {
        eca (mV) 
	ica 		(mA/cm2) 
	minf hinf 	(1)
	mtau (ms) htau 	(ms) 
	alphah (/ms) betah	(/ms)
}
 
STATE {
	m h
}

BREAKPOINT { 
	SOLVE states METHOD cnexp 
	ica = gbar * m * m * h * ( v - eca ) 
	alphah = hinf/htau
	betah = 1/htau - alphah
}
 
INITIAL { 
	settables(v) 

	h  = hinf
	m  = 0
} 

DERIVATIVE states { 
	settables(v) 
	m' = ( minf - m ) / mtau 
	h' = ( hinf - h ) / htau
}

UNITSOFF 

PROCEDURE settables(v(mV)) { 
	TABLE minf, mtau,hinf, htau FROM -120 TO 40 WITH 641
        minf  = 1 / ( 1 + exp( ( -v - 52 ) / 7.4 ) )
        mtau  = 1 + .33 / ( exp( ( v + 27.0 ) / 10.0 ) + exp( ( - v - 102 ) / 15.0 ) )

        hinf  = 1 / ( 1 + exp( ( v + 80 ) / 5 ) )
        htau = 28.30 +.33 / (exp(( v + 48.0)/ 4.0) + exp( ( -v - 407.0) / 50.0 ) )

}

UNITSON
UNITS { 
	(mV) = (millivolt) 
	(mA) = (milliamp) 
} 

NEURON { 
	SUFFIX IKs
	USEION k READ ek WRITE ik
	RANGE gbar, ik, tauM, kh
}

PARAMETER { 
	gbar =  0.0 	(mho/cm2)
	ek   = -70	(mV)
	kh    =  1.5       
	sha  =  0   (mV)
	shi  =  -3   (mV) 
	tauM = 10   (ms) 
} 

ASSIGNED { 
    v  (mV)
	ik  		(mA/cm2) 
	minf 		
	mtau 		(ms)
	hinf
	htau		(ms)
} 

STATE {
	m h
}

INITIAL { 
	Rates(v) 
	m = minf
	h = hinf
} 

BREAKPOINT { 
	SOLVE states METHOD cnexp
	ik = gbar * m * h * ( v - ek )
} 


DERIVATIVE states { 
	Rates(v) 
	m' = ( minf - m ) / mtau
	h' = ( hinf - h ) / htau 
}

UNITSOFF
 
PROCEDURE Rates(v) { 
	
	minf  = 1 / ( 1 + exp( -(v + 34 - sha) / 6.5 ) ) 
	mtau  = tauM 
	hinf = 1 / ( 1 + exp((v + 65 - shi) / 6.6 ))
	htau = 200 + kh*220/( 1 + exp( -(v + 71.6) / 6.85 ))
}
UNITSON
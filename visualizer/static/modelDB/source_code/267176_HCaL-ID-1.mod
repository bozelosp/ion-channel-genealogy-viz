TITLE Calcium high-threshold L type current for RD Traub, J Neurophysiol 89:909-921, 2003

COMMENT

	Implemented by Maciej Lazarewicz 2003 (mlazarew@seas.upenn.edu)
ENDCOMMENT

INDEPENDENT { t FROM 0 TO 1 WITH 1 (ms) }

UNITS { 
	(mV) = (millivolt) 
	(mA) = (milliamp) 
}
 
NEURON { 
	SUFFIX HCaL
	USEION ca WRITE ica
	RANGE  gbar, ica
	RANGE alpha, beta : added to write rates for comparison with F
}

PARAMETER { 
	gbar = 0.0 	(mho/cm2)
	v  		(mV)  
}
 
ASSIGNED { 
	ica 		(mA/cm2) 
	alpha (/ms)
  beta	(/ms)
}
 
STATE {
	m
}

BREAKPOINT { 
	SOLVE states METHOD cnexp
	ica = gbar * m^3* ( v - 80 )  
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

PROCEDURE settables(v (mV)) { LOCAL tmp
	TABLE alpha, beta FROM -120 TO 40 WITH 641

	 alpha = 0.025/(1+exp((v+30+10)/-7))
   

   beta=0.025/(1+exp((v+30+10)/7))
}

UNITSON

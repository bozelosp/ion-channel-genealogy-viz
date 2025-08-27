VERBATIM
#include <stdlib.h> 
ENDVERBATIM
 
UNITS {
	(mA) =(milliamp)
	(mV) =(millivolt)
	(uF) = (microfarad)
	(molar) = (1/liter)
	(nA) = (nanoamp)
	(mM) = (millimolar)
	(um) = (micron)
	FARADAY = 96520 (coul)
	R = 8.3134	(joule/degC)
}
 
NEURON { 
	SUFFIX ch_HCN 
	NONSPECIFIC_CURRENT i
	RANGE gmax, g, i, e
	RANGE hinf
	RANGE slow_tau 
	RANGE myi
}
 
 
PARAMETER {
	gmax  (mho/cm2)
	e (mV)
}
 
STATE {
	h
}
 
ASSIGNED {
	v (mV) 
	celsius (degC)
	dt (ms)    
  
	g (mho/cm2)
 	i (mA/cm2)
	
	hinf 
 	
 	slow_tau (ms) 
 	myi (mA/cm2)
} 

BREAKPOINT {
	SOLVE states METHOD cnexp
		
	g = gmax*h*h
	i = g*(v - e)
	myi = i
}
 
UNITSOFF
 
INITIAL { 
	trates(v)
	h = hinf
}

DERIVATIVE states {	
	trates(v)
	h' = (hinf-h)/slow_tau 
}
 
LOCAL q10
PROCEDURE trates(v) {  
	TABLE hinf, slow_tau 
	DEPEND celsius 
	FROM -120 TO 100 WITH 220
                           
    
    q10 = 3^((celsius - 34)/10)
       
	hinf =  1 / (1 + exp( (v+91)/10 ))

	
	

	
	slow_tau = (80*1.5 + .75*172.7 / (1+exp(-(v+59.3)/-0.83)))/q10

}

UNITSON
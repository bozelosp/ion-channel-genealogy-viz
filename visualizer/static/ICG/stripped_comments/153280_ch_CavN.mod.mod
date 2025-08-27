VERBATIM
#include <stdlib.h> 
ENDVERBATIM
 
UNITS {
	(mA) =(milliamp)
	(mV) =(millivolt)
	(molar) = (1/liter)
	(mM) = (millimolar)
	FARADAY = 96520 (coul)
	R = 8.3134	(joule/degC)
}
 
 
NEURON {
	SUFFIX ch_CavN				
	USEION ca READ eca WRITE ica VALENCE 2 
	RANGE g
	RANGE gmax
	RANGE cinf, ctau, dinf, dtau
	RANGE myi
	THREADSAFE
}
 
INDEPENDENT {t FROM 0 TO 100 WITH 100 (ms)}


PARAMETER {
	v (mV) 					
      celsius (degC) 
	gmax =1.0 (mho/cm2)		
}
 
STATE {
	c d		
}
 
ASSIGNED {			
	dt (ms) 				

	ica (mA/cm2)	
	g (mho/cm2)	
	eca (mV)		

	cinf dinf
	ctau (ms)
	dtau (ms) 
	cexp dexp      
	myi (mA/cm2)
}

BREAKPOINT {
	SOLVE states 
    g = gmax*c*c*d
	ica = g*(v-eca)
	myi = ica
}
 
UNITSOFF
 
INITIAL {
	trates(v)
	c = cinf
	d = dinf
}

? states 
PROCEDURE states() {	
        trates(v)	
	c = c + cexp*(cinf-c)
	d = d + dexp*(dinf-d)
        
}
 
LOCAL q10

PROCEDURE rates(v) {  
                      
        LOCAL  alpha, beta, sum
       
       q10 = 3^((celsius - 34)/10)
                
        alpha = -0.19*vtrap(v-19.88,-10)
	beta = 0.046*exp(-v/20.73)
	sum = alpha+beta        
	ctau = 1/sum      cinf = alpha/sum
                
	alpha = 0.00016*exp(-v/48.4) 
	beta = 1/(exp((-v+39)/10)+1)
	sum = alpha+beta        
	dtau = 1/sum      dinf = alpha/sum
}

PROCEDURE trates(v) {  
                      
	LOCAL tinc
	TABLE  cinf, cexp, dinf, dexp, ctau, dtau
	DEPEND dt, celsius FROM -100 TO 100 WITH 200
                           
	rates(v)	
				
				

	tinc = -dt * q10
	cexp = 1 - exp(tinc/ctau)
	dexp = 1 - exp(tinc/dtau)
}

FUNCTION vtrap(x,y) {  
	if (fabs(x/y) < 1e-6) {
		vtrap = y*(1 - x/y/2)
	}else{  
		vtrap = x/(exp(x/y) - 1)
	}
}

UNITSON
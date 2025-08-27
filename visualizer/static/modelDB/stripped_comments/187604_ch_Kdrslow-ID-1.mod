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
	SUFFIX ch_Kdrslow 
	USEION k READ ek WRITE ik  VALENCE 1
	RANGE g, gmax, ninf, ntau, ik
	RANGE myi
	THREADSAFE
}
 
PARAMETER {
	v (mV) 
	celsius (degC) 
	dt (ms) 

	ek  (mV)
	gmax (mho/cm2)
}
 
STATE {
	n
}
 
ASSIGNED {		     
	g (mho/cm2)
	ik (mA/cm2)
	ninf
	ntau (ms)
	nexp
	myi (mA/cm2)
} 

BREAKPOINT {
	SOLVE states
	g = gmax*n*n*n*n
	ik = g*(v-ek)
	myi = ik
}

UNITSOFF
 
INITIAL {
	trates(v)

	n = ninf
}

PROCEDURE states() {	
	trates(v)	
	n = n + nexp*(ninf-n)
}
 
LOCAL q10
PROCEDURE rates(v) {  
                      
	LOCAL  alpha, beta, sum
	q10 = 3^((celsius - 34)/10)
	

	
	alpha = -0.028*vtrap((v+65-35),-6)
	beta = 0.1056/exp((v+65-10)/40)
	sum = alpha+beta        
	ntau = 1/sum
	ninf = alpha/sum
	
}
 
PROCEDURE trates(v) {  
                      
	LOCAL tinc
	TABLE  ninf, nexp, ntau
	DEPEND dt, celsius
	FROM -100 TO 100 WITH 200
							   
	rates(v)	
	
	

	tinc = -dt * q10
	nexp = 1 - exp(tinc/ntau)
}
 
FUNCTION vtrap(x,y) {  
        if (fabs(x/y) < 1e-6) {
                vtrap = y*(1 - x/y/2)
        }else{  
                vtrap = x/(exp(x/y) - 1)
        }
}
 
UNITSON
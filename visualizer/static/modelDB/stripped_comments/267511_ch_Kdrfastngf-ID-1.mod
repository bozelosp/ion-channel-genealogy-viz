NEURON { 
	SUFFIX ch_Kdrfastngf
	USEION k READ ek WRITE ik VALENCE 1
	RANGE g, gmax, ninf, ntau, ik
	RANGE myi, offset5, offset6, slope5, slope6
	THREADSAFE
}

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
 
PARAMETER {

	
	gmax (mho/cm2)
	offset5=10 (mV)
	offset6=10 (mv)
	slope5=.07 (1)
	slope6=.264 (1)
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
	ek (mV)
	v (mV) 
	celsius (degC) 
	dt (ms) 
} 

BREAKPOINT {
	SOLVE states
	g = gmax*n*n*n*n
	ik = g*(v-ek)
	myi =  ik
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
                      
	LOCAL  alpha, beta, sum, tinc
	
	q10 = 3^((celsius - 34)/10)

	
	alpha = -1*slope5*vtrap((v+65-47-offset5),-6)
	beta = slope6/exp((v+65-22-offset6)/40)
	sum = alpha+beta        
	ntau = 1/sum
	ninf = alpha/sum	
	
	tinc = -dt * q10
	nexp = 1 - exp(tinc/ntau)
}
 
PROCEDURE trates(v) {  
                      
	LOCAL tinc
	TABLE ninf, nexp, ntau
	DEPEND dt, celsius, slope5, slope6, offset5, offset6
	FROM -100 TO 100 WITH 200
						   
	rates(v)	
	
	

	
	
}
 
FUNCTION vtrap(x,y) {  
        if (fabs(x/y) < 1e-6) {
                vtrap = y*(1 - x/y/2)
        }else{  
                vtrap = x/(exp(x/y) - 1)
        }
}
 
UNITSON
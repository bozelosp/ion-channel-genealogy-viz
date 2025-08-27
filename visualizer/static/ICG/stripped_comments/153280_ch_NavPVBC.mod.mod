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
	SUFFIX ch_NavPVBC 
	USEION na READ ena WRITE ina VALENCE 1
	RANGE g, gmax, minf, mtau, hinf, htau, ina, m, h
	RANGE myi, vshift
	THREADSAFE
}
 
PARAMETER {
	ena  (mV)
	gmax (mho/cm2)  
	vshift (mV)
}
 
STATE {
	m h
}
 
ASSIGNED {
	v (mV) 
	celsius (degC) 
	dt (ms) 

	g (mho/cm2)
	ina (mA/cm2)
	minf
	hinf
	mtau (ms)
	htau (ms)
	mexp
	hexp 
	myi (mA/cm2)
} 

BREAKPOINT {
	SOLVE states
	g = gmax*m*m*m*h  
	ina = g*(v - ena)
	myi = ina
}
 
UNITSOFF
 
INITIAL {
	trates(v)
	m = minf
	h = hinf
}

PROCEDURE states() {	
	trates(v)			
	m = m + mexp*(minf-m)
	h = h + hexp*(hinf-h)
}
 
LOCAL q10	
PROCEDURE rates(v) {  
                      
	LOCAL  alpha, beta, sum	

	q10 = 3^((celsius - 34)/10) 

	
	alpha = -0.3*vtrap((v+60-17+vshift),-5)
	beta = 0.3*vtrap((v+60-45+vshift),5)
	sum = alpha+beta        
	mtau = 1/sum 
	minf = alpha/sum
	
	
	alpha = 0.23/exp((v+60+5+vshift)/20)
	beta = 3.33/(1+exp((v+60-47.5+vshift)/-10))
	sum = alpha+beta
	htau = 1/sum 
	hinf = alpha/sum 	
}
 
PROCEDURE trates(v) {  
                      
	LOCAL tinc	
	TABLE minf, mexp, hinf, hexp, mtau, htau
	DEPEND dt, celsius FROM -100 TO 100 WITH 200
                                   
	rates(v)	
				
				

	tinc = -dt * q10

	mexp = 1 - exp(tinc/mtau)
	hexp = 1 - exp(tinc/htau)
 }
 
FUNCTION vtrap(x,y) {  
	if (fabs(x/y) < 1e-6) {
		vtrap = y*(1 - x/y/2)
	}else{  
		vtrap = x/(exp(x/y) - 1)
	}
}
 
UNITSON
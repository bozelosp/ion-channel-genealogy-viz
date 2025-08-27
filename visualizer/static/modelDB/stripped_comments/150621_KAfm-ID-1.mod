UNITS {
        (mA) = (milliamp)
        (mV) = (millivolt)
}
 
NEURON {
 	SUFFIX KAfm
	USEION k WRITE ik
	RANGE gkafmbar, gkafm
}
 
INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}
 
PARAMETER {
  	
	ek	= -73	(mV)
	gkafmbar= 0.00009 (mho/cm2) 
	Etemp	= 22 
	Vsm	= -33.1
	ksm	= 7.5
	tom	= 1
	Vsh	= -70.4
	ksh	= -7.6
	toh	= 25	     
}
 
STATE {
        m h
}
 
ASSIGNED {
        ik (mA/cm2)
 	v  (mV)
	celsius		(degC)
	minf
	hinf
	mtau
	htau
	gkafm
	
}
 
BREAKPOINT {
        SOLVE states METHOD cnexp
        gkafm = gkafmbar*m*h
        ik = gkafm*(v - ek)
  
}
 
UNITSOFF
 
INITIAL {
	rates(v)
	m= minf
	h= hinf
}

DERIVATIVE states {  
        rates(v)      
       
	m' = ( minf - m ) / mtau
	h' = (hinf - h ) / htau
}
 
PROCEDURE rates(v) {  
        LOCAL  q10,tadj  
        q10 = 2.5
	tadj=q10^((celsius-Etemp)/10)
        minf=1/(1+exp(-(v-Vsm)/ksm))
	hinf=1/(1+exp(-(v-Vsh)/ksh))
	mtau=tom/tadj
	htau=toh/tadj	      
}
 
UNITSON
UNITS {
        (mA) = (milliamp)
        (mV) = (millivolt)
}
 
NEURON {
 	SUFFIX Kirm
	USEION k WRITE ik
	RANGE gkirmbar, gkirm, minf
}
 
INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}
 
PARAMETER {
  	
	ek	= -90	(mV)
	gkirmbar= 0.00015 (mho/cm2) 
	tau	= 0.01	 (ms)
	Vsm	= -100
	ksm	= -10   
}
 
STATE {
        m
}
 
ASSIGNED {
	v  (mV)
        ik (mA/cm2)
	celsius		(degC)
 	minf
	mtau
	gkirm
}
 
BREAKPOINT {
        SOLVE states METHOD cnexp
        gkirm = gkirmbar*m
        ik = gkirm*(v - ek)
  
}
 
UNITSOFF
 
INITIAL {
	rates(v)
	m= minf
}

DERIVATIVE states {  
        rates(v)      
       
	m' = ( minf - m ) / mtau

}
 
PROCEDURE rates(v) {  
                      
        
        minf=1/(1+exp(-(v-Vsm)/ksm))
	mtau=tau     
}
 
UNITSON
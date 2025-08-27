UNITS {
        (mA) = (milliamp)
        (mV) = (millivolt)
}
 
NEURON {
 	SUFFIX NaPm
	USEION na WRITE ina
	RANGE gnapmbar, gnapm, minf, mtau
}
 
INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}
 
PARAMETER {
  	
	ena	= 45	(mV)
	gnapmbar= 0.00002 (mho/cm2) 
	Etemp	= 22	
	Vsm	= -47.8
	ksm	= 3.1
	tom	= 1
	
		     
}
 
STATE {
        m
}
 
ASSIGNED {
	v  (mV)
        ina (mA/cm2)
	celsius		(degC)
 	minf
	mtau
	gnapm
}
 
BREAKPOINT {
        SOLVE states METHOD cnexp
        gnapm = gnapmbar*m
        ina = gnapm*(v - ena)
  
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
        LOCAL  q10, tadj
        q10 = 2.5
	tadj=q10^((celsius-Etemp)/10)
        minf=1/(1+exp(-(v-Vsm)/ksm))
	mtau=tom/tadj
      
}
 
UNITSON
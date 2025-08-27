UNITS {
        (mA) = (milliamp)
        (mV) = (millivolt)
}
 
NEURON {
 	SUFFIX KAsm
	USEION k READ ek WRITE ik
	RANGE gkasmbar, gkasm, minf, hinf, mtau, htau
}
 
INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}
 
PARAMETER {
  	
	
	gkasmbar= 0.00032 (mho/cm2) 
	Etemp	= 22  

	Vsm	= -25.6
	ksm	= 13.3
	Vsh	= -78.8
	ksh	= -10.4
	tom	= 131.4
	Vtm	= -37.4
	ktm	= 27.3
	Vth	= -38.2
	kth	= 28
	hint	= 0.46     

}
 
STATE {
        m h
}
 
ASSIGNED {
        ek (mV)
	v  (mV)
        ik (mA/cm2)
	celsius		(degC)
 	minf
	hinf
	mtau
	htau
	gkasm
}
 
BREAKPOINT {
        SOLVE states METHOD cnexp
        gkasm = gkasmbar*m*h
        ik = gkasm*(v - ek)
  
}
 
UNITSOFF
 
INITIAL {
	rates(v)
	m= minf
	h= hint
}

DERIVATIVE states {  
        rates(v)      
       
	m' = ( minf - m ) / mtau
	h' = (hinf - h ) / htau
}
 
PROCEDURE rates(v) {  
        LOCAL  q10, tadj
        q10 = 2.5
	tadj=q10^((celsius-Etemp)/10)
        minf=1/(1+exp(-(v-Vsm)/ksm))
	hinf=1/(1+exp(-(v-Vsh)/ksh))
	mtau=tom/(exp(-(v-Vtm)/ktm)+exp((v-Vtm)/ktm))/tadj
	htau=(1790+2930*exp(-((v-Vth)/kth)^2)*((v-Vth)/kth))/tadj
	      
}
 
UNITSON
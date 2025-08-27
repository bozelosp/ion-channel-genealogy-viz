NEURON {
	SUFFIX car
	USEION ca READ cai, eca WRITE ica 
        RANGE gcabar, ica, po
	GLOBAL hinf, minf
}

UNITS {
	(mA) = (milliamp)
	(mV) = (millivolt)
	(molar) = (1/liter)
	(mM) =	(millimolar)
}

PARAMETER {           

	gcabar = 1.0   (mho/cm2)  
	zetam = -3.4		
	zetah = 2		
	vhalfm =-21 (mV)	
	vhalfh =-40 (mV)
	tm0=1.5(ms)
	th0=75(ms)
}



ASSIGNED {     
	v            (mV)
	celsius      (degC)
	ica          (mA/cm2)
	po
	cai          (mM)       
	eca             (mV)
        minf
        hinf
}



STATE {	
	m 
	h 
}  

INITIAL {
	rates(v)
        m = minf
        h = hinf
}

BREAKPOINT {
	SOLVE states METHOD cnexp
	po = m*m*h
	ica = gcabar *po* (v-eca)

}


DERIVATIVE states {
	rates(v)
	m' = (minf -m)/tm0
	h'=  (hinf - h)/th0
}


PROCEDURE rates(v (mV)) { 
        LOCAL a, b
        
	a = alpm(v)
	minf = 1/(1+a)
        
        b = alph(v)
	hinf = 1/(1+b)
}



FUNCTION alpm(v(mV)) {
UNITSOFF
  alpm = exp(1.e-3*zetam*(v-vhalfm)*9.648e4/(8.315*(273.16+celsius))) 
UNITSON
}

FUNCTION alph(v(mV)) {
UNITSOFF
  alph = exp(1.e-3*zetah*(v-vhalfh)*9.648e4/(8.315*(273.16+celsius))) 
UNITSON
}
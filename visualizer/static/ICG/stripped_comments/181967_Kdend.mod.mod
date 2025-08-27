UNITS {
        (mA) = (milliamp)
        (mV) = (millivolt)
}
 
NEURON {
        SUFFIX Kdend
        USEION k READ ek WRITE ik
        RANGE gkdend, ik
        GLOBAL ninf, nexp, ntau
}
 
INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}
 
PARAMETER {
        v (mV)
        celsius = 24 (degC)
        dt (ms)
        gkdend = .0230 (mho/cm2)
        
}
 
STATE {
        n 
}
 
ASSIGNED {
        ik (mA/cm2)
        ninf 
	nexp 
	ntau (ms)
	ek (mV)
}
 
INITIAL {
	n = ninf
}

BREAKPOINT {
        SOLVE states
	ik = gkdend*n*n*n*n*(v - ek)    
}

PROCEDURE states() {	
	evaluate_fct(v)
	n = n + nexp*(ninf - n)
	VERBATIM
	return 0;
	ENDVERBATIM 
}
UNITSOFF
PROCEDURE evaluate_fct(v(mV)) {  
		      
                      
        LOCAL q10, tinc, alpha, beta
        TABLE ninf, nexp, ntau DEPEND dt, celsius FROM -200 TO 
100 WITH 300

		q10 = 1	
		tinc = -dt*q10
		alpha = 0.018*vtrap(-(v-20),21)
		beta = 0.0036*vtrap(v-30,12)
		ntau = 1/(alpha + beta)
		ninf = alpha*ntau
		nexp = 1-Exp(tinc/ntau)
}
FUNCTION vtrap(x,y) {	
		if (fabs(x/y) < 1e-6) {
			vtrap = y*(1 - x/y/2)
		}else{
			vtrap = x/(Exp(x/y) - 1)
		}
}
FUNCTION Exp(x) {
		if (x < -100) {
			Exp = 0
		}else{
			Exp = exp(x)
		}
}
UNITSON
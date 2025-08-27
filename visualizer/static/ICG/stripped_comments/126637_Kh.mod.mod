UNITS {
        (mA) = (milliamp)
        (mV) = (millivolt)
}
 
NEURON {
        SUFFIX Kh
	USEION k READ ek WRITE ik
        RANGE  gkbar, gk, minf, mexp, nexp, ik
} 
 
INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}
 
PARAMETER {
        v (mV)
        celsius = 37 (degC)
        dt (ms)
        gkbar	= .0003 (mho/cm2)
        
	mon = 1
	hon = 1
}
 
STATE {
        m 
}
 
ASSIGNED {
        ek (mV)
        ik (mA/cm2)
        gk minf  mexp nexp
}
 
BREAKPOINT {
        SOLVE states
        gk = gkbar *m
	ik = gk* (v-ek)
}
 
UNITSOFF
 
INITIAL {
	rates(v)
	m = minf
}

PROCEDURE states() {  
        rates(v)      
        m = mon * (0.8*(m + mexp*(minf-m)))
		+ hon * (0.2*(m + nexp*(minf-m)))
}
 
PROCEDURE rates(v) {  
                      
        LOCAL  q10, tinc
        TABLE minf, mexp,nexp DEPEND dt, celsius FROM -400 TO 300 WITH 700
        q10 = 3^((celsius - 37)/10)
        tinc = -dt * q10
                
        minf = 1/(1+exp((v+78)/7))
        mexp = 1 - exp(tinc/38)
        nexp = 1 - exp(tinc/319)
}

 
UNITSON
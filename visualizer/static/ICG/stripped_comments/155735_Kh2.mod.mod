UNITS {
        (mA) = (milliamp)
        (mV) = (millivolt)
}
 
NEURON {
        SUFFIX Kh2
        USEION k READ ek WRITE ik
        RANGE  gkbar, gk, minf, nexp, ik
	GLOBAL ek
} 
 
INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}
 
PARAMETER {
        v (mV)
        celsius = 37 (degC)
        dt (ms)
        gkbar   = .0003 (mho/cm2)
        ek     (mV)

}
 
STATE {
        m 
}
 
ASSIGNED {
        ik (mA/cm2)
        gk minf nexp
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
        m = m + nexp*(minf-m)
}

PROCEDURE rates(v) {  
                      
        LOCAL  q10, tinc
        TABLE minf,nexp DEPEND dt, celsius FROM -100 TO 100 WITH 200
        q10 = 3^((celsius - 37)/10)
        tinc = -dt * q10
                
        minf = 0.2/(1+exp((v+82)/7))
        nexp = 1 - exp(tinc/36.8)
}

 
UNITSON
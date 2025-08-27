UNITS {
        (mA) = (milliamp)
        (mV) = (millivolt)
}
 
NEURON {
        SUFFIX KMnew2
        USEION k READ ek WRITE ik
        RANGE  gkbar, gk, minf,  mexp, ik
	GLOBAL ek
} 

INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}
 
PARAMETER {
        v (mV)
        celsius = 37 (degC)
        dt (ms)
        gkbar   = .00004 (mho/cm2)
        ek     (mV)

}
 
STATE {
        m 
}
 
ASSIGNED {
        ik (mA/cm2)
        gk minf  mexp 
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
        m = m + mexp*(minf-m)
}
 
PROCEDURE rates(v) {  
                      
        LOCAL  q10, tinc, sum
        TABLE minf, mexp DEPEND dt, celsius FROM -100 TO 100 WITH 200
        q10 = 3^((celsius - 37)/10)
        tinc = -dt * q10
                



        sum = (3.3*exp((v + 35)/20) + exp(-(v+35)/20))/200
        minf = 1.0 / (1+exp(-(v+35)/10))
        mexp = 1 - exp(tinc*sum)
               
}

 
UNITSON
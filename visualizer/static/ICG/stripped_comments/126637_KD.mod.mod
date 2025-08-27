UNITS {
        (mA) = (milliamp)
        (mV) = (millivolt)
}
 
NEURON {
        SUFFIX KD
        USEION k READ ek WRITE ik
        RANGE  gkbar, ik, gk, minf, hinf, mexp, hexp
} 
 
INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}
 
PARAMETER {
        v (mV)
        celsius = 37 (degC)
        dt (ms)
        gkbar = .0045 (mho/cm2)
        
	mon = 1
	hon = 1
}
 
STATE {
        m h
}
 
ASSIGNED {
        ek (mV)
        ik (mA/cm2)
        gk minf hinf mexp hexp 
}
 
BREAKPOINT {
        SOLVE states
        gk = gkbar * m*h
	ik = gk* (v-ek)
}
 
UNITSOFF
 
INITIAL {
	rates(v)
	m = minf
	h = hinf
}

PROCEDURE states() {  
        rates(v)      
        m = mon * (m + mexp*(minf-m))
        h = hon * (h + hexp*(hinf-h))
}
 
PROCEDURE rates(v) {  
                      
        LOCAL  q10, tinc, alpha, beta, sum
        TABLE minf, mexp, hinf, hexp DEPEND dt, celsius FROM -400 TO 300 WITH 700
        q10 = 3^((celsius - 37)/10)
        tinc = -dt * q10
                
        alpha = 8.5/(1+exp((v+17)/(-12.5)))
        beta =  35/(1+exp((v+99)/14.5))
        sum = alpha + beta
        minf = alpha/sum
        mexp = 1 - exp(tinc*sum/10)
                
        alpha = 0.0015/(1+exp((v+89)/8))
        beta = 0.0055/(1+exp((v+83)/(-8)))
        sum = alpha + beta
        hinf = alpha/sum
        hexp = 1 - exp(tinc*sum*1.6)
}

 
UNITSON
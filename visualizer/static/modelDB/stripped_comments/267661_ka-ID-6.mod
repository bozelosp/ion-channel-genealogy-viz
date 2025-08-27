UNITS {
        (mA) = (milliamp)
        (mV) = (millivolt)
}
 
NEURON {
        SUFFIX ka
	USEION k WRITE ik
        RANGE  gkbar, gk, minf, hinf, mexp, hexp, ik, vshift_m, vshift_h, localtemp
} 
 
INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}
 
PARAMETER {
        v (mV)
        localtemp = 37
        dt (ms)
        gkbar	= .0012 (mho/cm2)
        ek	= -85 (mV)
		vshift_m=0
		vshift_h=0
}
 
STATE {
        m h
}
 
ASSIGNED {
        ik (mA/cm2)
        gk minf hinf mexp hexp
        q10 
}
 
BREAKPOINT {
        SOLVE states
        gk = gkbar *m*m*m*m*h 
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
        m = m + mexp*(minf-m)
        h = h + 0.04*hexp*(hinf-h)
}
 
PROCEDURE rates(v) {  
                      
        LOCAL  q10, tinc, alpha, beta, sum
        TABLE minf, mexp, hinf, hexp DEPEND dt, celsius FROM -100 TO 100 WITH 200
        q10 = 3^((localtemp - 37)/10)
        tinc = -dt * q10
                
        alpha = 1.4/(1+exp((v+27 + (-20))/(-12)))
        beta =  0.49/(1+exp((v+30+(-20))/4))
        sum = alpha + beta
        minf = alpha/sum
        mexp = 1 - exp(tinc*sum)
                
        alpha = 0.0175/(1+exp((v+50+(-20))/8))
        beta = 1.3/(1+exp((v+13+(-20))/(-10)))
        sum = alpha + beta
        hinf = alpha/sum
        hexp = 1 - exp(tinc*sum)
}

 
UNITSON
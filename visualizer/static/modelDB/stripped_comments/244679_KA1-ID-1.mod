UNITS {
        (mA) = (milliamp)
        (mV) = (millivolt)
}
 
NEURON {
        SUFFIX KA1
	USEION k WRITE ik
        RANGE  gkbar, gk, minf, hinf, mexp, hexp, ik, vh1, vm1, vm2, t1,t2
} 
 
INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}
 
PARAMETER {
        v (mV)
        celsius = 37 (degC)
        dt (ms)
        gkbar	= .015 (mho/cm2)
        ek	= -85 (mV)
	mon = 1
	hon = 1
    vh1=85 (mV)
    vm1 (mV)
    vm2 (mV)
    t1=0.3
    t2=50
}
 
STATE {
        m h
}
 
ASSIGNED {
        ik (mA/cm2)
        gk minf hinf mexp hexp 
}
 
BREAKPOINT {
        SOLVE states
        gk = gkbar *m*m*m* m*h 
	ik = gk* (v-ek)
}
 
UNITSOFF
 
INITIAL {
	rates(v,vm1,vm2)
	m = minf
	h = hinf
}

PROCEDURE states() {  
        rates(v,vm1,vm2)      
        m = mon * (m + mexp*(minf-m))
        h = hon * (h + hexp*(hinf-h))
}
 
PROCEDURE rates(v, vm1,vm2) {  
                      
        LOCAL  q10, tinc, alpha, beta, sum
        
        q10 = 3^((celsius - 37)/10)
        tinc = -dt * q10
                
        alpha = 1.4/(1+exp((v+vm1)/(-12)))
        beta =  0.49/(1+exp((v+vm2)/4))
        sum = alpha + beta
        minf = (alpha/sum)*(1/(1+exp(-t1*(v+t2))))
        mexp = 1 - exp(tinc*sum)
                
        alpha = 0.0175/(1+exp((v+vh1)/8))
        beta = 1.3/(1+exp((v+13)/(-10)))
        sum = (alpha + beta)
        hinf = alpha/sum
       
        hexp = 1 - exp(tinc*sum)
}

 
UNITSON
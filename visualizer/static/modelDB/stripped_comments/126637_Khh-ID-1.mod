UNITS {
        (mA) = (milliamp)
        (mV) = (millivolt)
}
 
NEURON {
        SUFFIX Khh
        USEION k WRITE ik
        RANGE   gk,  gkbar, ik
        GLOBAL  ninf, nexp
}
 
INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}
 
PARAMETER {
        v (mV)
        celsius = 37 (degC)
        dt (ms)
        gkbar = .036 (mho/cm2)
        ek = -85(mV)
	non = 1
}
 
STATE {
         n
}
 
ASSIGNED {
        ik (mA/cm2)
        gk ninf nexp
}
 
BREAKPOINT {
        SOLVE states
        gk  = gkbar*n*n*n*n

        ik = gk*(v - ek)      
}
 
UNITSOFF
 
INITIAL {
	rates(v)
	n = ninf
}

PROCEDURE states() {  
        rates(v)      
        n = non * (n + nexp*(ninf-n))
}
 
PROCEDURE rates(v) {  
                      
        LOCAL  q10, tinc, alpha, beta, sum
        TABLE ninf, nexp DEPEND dt, celsius FROM -400 TO 300 WITH 700
        q10 = 3^((celsius - 37)/10)
        tinc = -dt * q10
                
        alpha = .01*vtrap(-(v+55),10) 
        beta = .125*exp(-(v+65)/80)
        sum = alpha + beta
        ninf = alpha/sum
        nexp = 1 - exp(tinc*sum)
}

FUNCTION vtrap(x,y) {  
        if (fabs(x/y) < 1e-6) {
                vtrap = y*(1 - x/y/2)
        }else{
                vtrap = x/(exp(x/y) - 1)
        }
}
 
UNITSON
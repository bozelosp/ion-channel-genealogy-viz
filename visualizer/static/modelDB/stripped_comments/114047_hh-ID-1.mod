UNITS {
        (mA) = (milliamp)
        (mV) = (millivolt)
}
 
NEURON {
        SUFFIX HH
        USEION na READ ena WRITE ina
        USEION k READ ek WRITE ik
        RANGE gnabar, gkbar
        GLOBAL minf, hinf, ninf, mexp, hexp, nexp
}
 
PARAMETER {
        v (mV)
        celsius = 5 (degC)
        dt (ms)
        gnabar = .184 (mho/cm2)
        ena = 55 (mV)
        gkbar = .14 (mho/cm2)
        ek = -90 (mV)
}
 
STATE {
        m h n
}
 
ASSIGNED {
        ina (mA/cm2)
        ik (mA/cm2)
        minf hinf ninf mexp hexp nexp
}
 
BREAKPOINT {
        SOLVE states
        ina = gnabar*m*m*m*h*(v - ena)
        ik = gkbar*n*n*n*n*(v - ek)      
}
 
UNITSOFF
 
INITIAL {
	rates(v)
	m = minf
	h = hinf
	n = ninf
}

PROCEDURE states() {  
        rates(v)      
        m = m + mexp*(minf-m)
        h = h + hexp*(hinf-h)
        n = n + nexp*(ninf-n)
}
 
PROCEDURE rates(v) {  
                      
        LOCAL  q10, tinc, alpha, beta, sum
        TABLE minf, mexp, hinf, hexp, ninf, nexp DEPEND dt, celsius FROM -100 TO 100 WITH 200
        q10 = 3^((celsius - 6.3)/10)
        tinc = -dt * q10
                
        alpha = .1 * vtrap(-(v+35),10)
        beta =  4 * exp(-(v+60)/18)
        sum = alpha + beta
        minf = alpha/sum
        mexp = 1 - exp(tinc*sum)
                
        alpha = .07 * exp(-(v+58)/20)
        beta = 1 / (exp(-(v+28)/10) + 1)
        sum = alpha + beta
        hinf = alpha/sum
        hexp = 1 - exp(tinc*sum)
                
        alpha = .01*vtrap(-(v+34),10) 
        beta = .125*exp(-(v+44)/80)
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
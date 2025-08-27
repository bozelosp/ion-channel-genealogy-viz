UNITS {
        (mA) = (milliamp)
        (mV) = (millivolt)
}
 
NEURON {
        SUFFIX CA1CaG
        USEION ca WRITE ica
        RANGE gcabar,gca,ica
        GLOBAL minf, hinf, mexp, hexp, tc
}
 
INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}
 
PARAMETER {
        v (mV)
        celsius = 36 (degC)
        dt (ms)
        gcabar = 0.01 (mho/cm2)
        eca = 80 (mV)
        tc = 1 (1)
}
 
STATE {
        m h
}
 
ASSIGNED {
        ica (mA/cm2)
        minf hinf mexp hexp
}
 
BREAKPOINT {
        SOLVE states
        ica = gcabar*m*m*h*(v - eca)
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
        h = h + hexp*(hinf-h)
}
 
PROCEDURE rates(v) {  
                      
        LOCAL  alpha, beta, sum
        TABLE minf, mexp, hinf, hexp DEPEND dt, tc FROM -100 TO 100 WITH 200
    
        alpha = -0.16 * vtrap(v+26,-4.5)
        beta =  0.04 * vtrap(v+12,10)
        sum = alpha + beta
        minf = alpha/sum
        mexp = 1 - exp(-dt*sum)
    
        alpha = 2 / exp((v+94)/10)
        beta = 8 / (exp(-(v-68)/27) + 1)
        sum = alpha + beta
        hinf = alpha/sum
        hexp = 1 - exp(-dt*sum/tc)
}
 
FUNCTION vtrap(x,y) {  
        if (fabs(x/y) < 1e-6) {
                vtrap = y*(1 - x/y/2)
        }else{
                vtrap = x/(exp(x/y) - 1)
        }
}
 
UNITSON
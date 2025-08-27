UNITS {
        (mA) = (milliamp)
        (mV) = (millivolt)
         (S) = (siemens)
}
 
? interface
NEURON {
        SUFFIX hhin
        USEION na READ ena WRITE ina
        USEION k READ ek WRITE ik
        NONSPECIFIC_CURRENT il
        RANGE gnabar, gkbar, gl, el, gna, gk, a
        GLOBAL minf, hinf, ninf, sinf, mtau, htau, stau, ntau
}
 
PARAMETER {
        gnabar = 0.0   (S/cm2)    <0,1e9>
        gkbar  = .036  (S/cm2)    <0,1e9>
        gl     = 0.0 (S/cm2)    <0,1e9>
        el     = -54.  (mV)
        a      = 1     (1)
}
 
STATE {
        m h n s
}
 
ASSIGNED {
        v       (mV)
        celsius (degC)
        ena     (mV)
        ek      (mV)

        gna (S/cm2)
        gk  (S/cm2)
        ina (mA/cm2)
        ik  (mA/cm2)
        il  (mA/cm2)
        minf hinf sinf ninf
    mtau (ms) htau (ms) stau (ms) ntau (ms)
}
 

 
? currents
BREAKPOINT {
        SOLVE states METHOD cnexp
        gna = gnabar*m*m*m*h*s
        ina = gna*(v - ena)
        gk = gkbar*n*n*n*n
        ik = gk*(v - ek)      
        il = gl*(v - el)
}
 
 
INITIAL {
    rates(v)
    m = minf
    h = hinf
    s = sinf
    n = ninf
}

? states
DERIVATIVE states {  
        rates(v)
        m' =  (minf-m)/mtau
        h' =  (hinf-h)/htau
        s' =  (sinf-s)/stau
        n' =  (ninf-n)/ntau
}
 
LOCAL q10


? rates
PROCEDURE rates(v(mV)) {  
                      
        LOCAL  alpha, beta, sum
        TABLE minf, mtau, hinf, sinf, htau, ninf, stau, ntau DEPEND celsius FROM -100 TO 100 WITH 200


UNITSOFF
        q10 = 3^((celsius - 6.3)/10)
                
        alpha = .1 * vtrap(-(v+40),10) 
        beta =  4 * exp(-(v+65)/18)
        sum = alpha + beta
    mtau = 1/(q10*sum)
        minf = alpha/sum
                
        alpha = .07 * exp(-(v+65)/20)
        beta = 1 / (exp(-(v+35)/10) + 1)
        sum = alpha + beta
    htau = 1/(q10*sum)
        hinf = alpha/sum
                
        alpha = 0.0077 / (1. + exp( (47.+v)/9. ) )
        beta  = 0.0077 / (1. + exp( -(47.+v)/9.  ) )
        sum = alpha + beta
    stau = 1/(sum)
        sinf = (1-a) + a*alpha/sum
                
        alpha = .01*vtrap(-(v+55),10) 
        beta = .125*exp(-(v+65)/80)
    sum = alpha + beta
        ntau = 1/(q10*sum)
        ninf = alpha/sum
}
 
FUNCTION vtrap(x,y) {  
        if (fabs(x/y) < 1e-6) {
                vtrap = y*(1 - x/y/2)
        }else{
                vtrap = x/(exp(x/y) - 1)
        }
}
 
UNITSON
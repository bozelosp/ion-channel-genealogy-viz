UNITS {
        (mA)    = (milliamp)
        (mV)    = (millivolt)
        (molar) = (1/liter)
        (mM)    = (millimolar)
        (S)     = (siemens)
        FARADAY = 96520 (coul)
        R       = (k-mole) (joule/degC)     
        KTOMV   = .0853 (mV/degC)
}



NEURON {
        SUFFIX iCaT
        USEION ca READ cai,cao WRITE ica VALENCE 2
        RANGE gcatbar, minf, mtau, hinf, htau, i
}



PARAMETER {                

        gcatbar = 0     (S/cm2) 
        ki  = 0.001     (mM)  
        cai = 5.e-5     (mM)      
        cao = 2         (mM)      
        tfa = 1                   
        tfi = 0.68
        tBase = 23.5    (degC)
        eca = 140       (mV)      
}



ASSIGNED {                       

        v               (mV)
        celsius         (degC)

        gcat            (S/cm2) 
        ica             (mA/cm2)
        i               (mA/cm2)
        asdf            

        minf
        mtau            (ms)
        hinf
        htau            (ms)
        }



STATE {  m h }                          



BREAKPOINT {

        SOLVE states METHOD cnexp
        gcat = gcatbar*m*m*h*h2(cai)    
        i   = gcat*ghk(v,cai,cao)       
        ica = i
}



INITIAL {                               
        rates(v)
        m = minf
        h = hinf
        gcat = gcatbar*m*m*h*h2(cai)
}



? states
DERIVATIVE states {
    
    rates(v)
    m' = (minf-m)/mtau
    h' = (hinf-h)/htau

}



PROCEDURE rates(v (mV)) { 

        LOCAL a

        UNITSOFF 
        a = alpm(v)
        mtau = 1.0 / (tfa*(a+betm(v)))  
        minf = a /( a + betm(v) )       

        a = alph(v)
        htau = 1.0 / (tfi*(a+beth(v)))  
        hinf = a / ( a + beth( v ) )    

        UNITSON
}



UNITSOFF

FUNCTION h2(cai(mM)) {
        h2 = ki/(ki+cai)
}


FUNCTION ghk(v(mV), ci(mM), co(mM)) (mV) {

        LOCAL nu,f

        f   = KTF(celsius)/2
        nu  = v/f
        ghk = -f*(1. - (ci/co)*exp(nu))*efun(nu)
}


FUNCTION KTF(celsius (degC)) (mV) { 
        KTF = ((25./293.15)*(celsius + 273.15))
}


FUNCTION efun(z) {
        if (fabs(z) < 1e-4) {
                efun = 1 - z/2
        }else{
                efun = z/(exp(z) - 1)
        }
}


FUNCTION alpm(v(mV)) {
    TABLE FROM -150 TO 150 WITH 200
    alpm = 0.1967*(-1.0*v+19.88)/(exp((-1.0*v+19.88)/10.0)-1.0)
}

FUNCTION betm(v(mV)) {
    TABLE FROM -150 TO 150 WITH 200
    betm = 0.046*exp(-v/22.73)
}

FUNCTION alph(v(mV)) {
    TABLE FROM -150 TO 150 WITH 200
    alph = 1.6e-2*exp(-(v+57.0)/19.0)
}

FUNCTION beth(v(mV)) {
    TABLE FROM -150 TO 150 WITH 200
    beth = 1.0/(exp((-v+15)/10)+1.0)
}

UNITSON
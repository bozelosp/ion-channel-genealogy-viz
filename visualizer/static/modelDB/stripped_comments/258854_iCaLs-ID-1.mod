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
        SUFFIX iCaLs
        USEION ca READ cai,cao WRITE ica
        RANGE gcalbar, minf, mtau, i
}



PARAMETER {                

        gcalbar = 0     (S/cm2) 
        ki  = 0.001     (mM)  
        cai = 5.e-5     (mM)      
        cao = 2         (mM)      
        tfa = 5                   
        eca = 140       (mV)      
}



ASSIGNED {                       

        v               (mV)
        celsius         (degC)

        gcal            (S/cm2) 
        ica             (mA/cm2)
        i               (mA/cm2)

        minf
        mtau            (ms)
        }



STATE {  m  }                      



BREAKPOINT {

        SOLVE states METHOD cnexp
        gcal = gcalbar*m*h2(cai) 
        i   = gcal*ghk(v,cai,cao)
        ica = i
}



INITIAL {                        
        rates(v)
        m = minf
        gcal = gcalbar*m*h2(cai)
}



? states
DERIVATIVE states {
    
    rates(v)
    m' = (minf-m)/mtau

}



PROCEDURE rates(v (mV)) { 

        LOCAL a

        UNITSOFF 
        a = alpm(v)
        mtau = 1.0 / (tfa*(a+betm(v)))  
        minf = a / (a+betm(v))          

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
        alpm = - 0.055 * (-27.01+v )/(exp( -(-27.01+v)/3.8) - 1)
}


FUNCTION betm(v(mV)) {
        TABLE FROM -150 TO 150 WITH 200
        betm =0.94*exp(-(63.01+v)/17)
}

UNITSON
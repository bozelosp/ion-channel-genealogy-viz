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
        SUFFIX iCaLd
        USEION ca READ cai,cao WRITE ica
        RANGE gcalbar, minf, mtau, hinf, htau, i
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
        hinf
        htau            (ms)
        }



STATE { 
        m
        h 
}                      



BREAKPOINT {

        SOLVE states METHOD cnexp
        gcal = gcalbar*m*m*m*h
        i    = gcal*(v - eca)
        ica  = i
}



INITIAL {                        

        rates(v)
        m   = 0
        h   = 1
}



? states
DERIVATIVE states {
    
    rates(v)
    m' = (minf-m)/mtau
    h' = (hinf-h)/htau

}



PROCEDURE rates(v (mV)) { 

        UNITSOFF 

        minf = 1.0 / (1.0 + exp((v+37.0)/(-1.0)))         
        mtau = 3.6

        hinf = 1.0 / (1.0 + exp((v+41.0)/(0.5)))          
        htau = 29.0

        UNITSON
}
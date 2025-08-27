UNITS {
    (mA) = (milliamp)
    (mV) = (millivolt)
}
 
? interface
NEURON {
    THREADSAFE
    SUFFIX KI                                                                                         
    USEION k READ ek WRITE ik                                                            
    RANGE gkbar, gk, scale_a                                                                                  
    GLOBAL ninf, ntau
}
 
PARAMETER {
    gkbar = .004 (mho/cm2)	 <0,1e9>   
    scale_a = 1.0
}
 
STATE {
    n                                
}
 
ASSIGNED {                                   
    v (mV)                          
    celsius (degC)
    ek (mV)
    
    gk (mho/cm2)               
    ik (milliamp/cm2)
    ninf
    ntau (ms)
}
 
 
? currents
BREAKPOINT {                                      
    SOLVE states METHOD cnexp
    gk = gkbar*n*n*n*n
    ik = gk*(v - ek)      
}
 
 
INITIAL {
    rates(v)
    n = ninf
}

DERIVATIVE states {  
    rates(v)
    n' = (ninf-n)/ntau
}
 
? rates
PROCEDURE rates(v(mV)) {  
                                                    
    LOCAL  alpha, beta, sum

UNITSOFF
    alpha = scale_a*.01*vtrap(-(v+55),10) 
    beta = scale_a*.125*exp(-(v+65)/80)
    sum = alpha + beta
    ntau = 1/sum
    ninf = alpha/sum
}
 
FUNCTION vtrap(x,y) {            
    if (fabs(x/y) < 1e-6) {
        vtrap = y*(1 - x/y/2)
    } else {
        vtrap = x/(exp(x/y) - 1)
    }
}
 
UNITSON
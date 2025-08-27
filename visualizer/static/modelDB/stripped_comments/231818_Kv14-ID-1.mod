UNITS {
    (mA) = (milliamp)
    (mV) = (millivolt)
}
 
NEURON {
    THREADSAFE
    SUFFIX Kv14                                                                                           
    USEION k READ ek WRITE ik                                                            
    RANGE gkbar, gk, scale_a, scale_i, ik                                                                                  
    GLOBAL vshift
}
 
PARAMETER {
    gkbar = .0 (mho/cm2)	 <0,1e9>   
    scale_a = 1.0
    scale_i = 1.0
	vshift = 0 (mV)
}
 
STATE {
    n h                                
}
 
ASSIGNED {                                   
    v (mV)                          
    celsius (degC)
    ek (mV)
    
    gk (mho/cm2)               
    ik (milliamp/cm2)
    ninf
    ntau (ms)
    hinf
    htau (ms) 
}
 
 
BREAKPOINT {                                      
    SOLVE states METHOD cnexp
    gk = gkbar*n*n*n*n*h
    ik = gk*(v - ek)      
}
 
 
INITIAL {
    rates(v)
    n = ninf
    h = hinf
}


DERIVATIVE states {  
    rates(v)
    n' = (ninf-n)/ntau
    h' = (hinf-h)/htau
}
 
PROCEDURE rates(v(mV)) {  
                                                    
    LOCAL  alpha, beta, sum

UNITSOFF
    alpha = scale_a*.01*vtrap((-55+vshift-v),10) 
    beta = scale_a*.125*exp((-65+vshift-v)/80)
    sum = alpha + beta
    ntau = 1/sum
    ninf = alpha/sum
    
    alpha = scale_i*0.0000256077*exp((vshift-v)/45.4217)
    beta = scale_i*0.0330402/(exp((-45.6599+vshift-v)/2.30235) + 1) 
    sum = alpha + beta
    htau = 1/sum
    hinf = alpha/sum
}
 
FUNCTION vtrap(x,y) {            
    if (fabs(x/y) < 1e-6) {
        vtrap = y*(1 - x/y/2)
    } else {
        vtrap = x/(exp(x/y) - 1)
    }
}
 
UNITSON
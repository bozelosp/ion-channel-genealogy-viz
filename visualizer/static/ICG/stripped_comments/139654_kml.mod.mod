UNITS {
        (mA) = (milliamp)
        (mV) = (millivolt)
        (nA) = (nanoamp)
}

NEURON {
        SUFFIX kml
        USEION k READ ek WRITE ik
        RANGE gbar, g, bn, gn, tn, ik
        GLOBAL ninf, ntau
}

INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}

PARAMETER {
        v (mV)
        ek = -90 (mV)
        gbar = 1.0 (S/cm2) <0,1e9>
        bn = -20 (mV)
        gn = 10 (mV) 
        tn = 3 (ms)
}

STATE {
        n
}

ASSIGNED {
    ik (mA/cm2) 
    g (S/cm2)
    ninf
    ntau (ms)
    }

BREAKPOINT {
	SOLVE states METHOD cnexp 

	g = gbar*n
    ik = g*(v - ek)

}

INITIAL {
    rates(v)
    n = ninf
}


DERIVATIVE states {  
        rates(v)
        n' =  (ninf-n)/ntau
}


UNITSOFF

PROCEDURE rates(v) {  
                      
    
    
    ninf = (1 / (1 + exp((bn - v) / (gn/2))))
    
    ntau = (2*tn) / ( exp((v-bn)/(2*gn)) + exp((-v+bn)/(2*gn)) )
}

UNITSON
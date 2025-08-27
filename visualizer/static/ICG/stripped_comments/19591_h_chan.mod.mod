UNITS {
        (mA) = (milliamp)
        (mV) = (millivolt)
}
 
NEURON {
        SUFFIX h_chan
        RANGE gbar
        NONSPECIFIC_CURRENT i
        GLOBAL minf, 
               am, bm, cm, dm, taum_min
	GLOBAL eh
}
 
INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}
 
PARAMETER {
        v (mV)
        dt (ms)
        gbar = 1.648e-4 (mho/cm2)
        
        

        am = 0.001
        bm = -0.1818
        cm = 0.5
        dm = -75
        taum_min = 1e-3

}
 
STATE {
        m
}
 
ASSIGNED {
        eh (mV)
        i (mA/cm2)
        minf
}
 
LOCAL mexp
 
BREAKPOINT {
        SOLVE states
        
	i = gbar*m*(v - eh)
}
 
UNITSOFF
 
INITIAL {
	rates(v)
	m = minf
}

PROCEDURE states() {  
        rates(v)      
        m = m + mexp*(minf-m)
        VERBATIM
        return 0;
        ENDVERBATIM
}
 
PROCEDURE rates(v) {  
                      
        LOCAL  tau,alpha,beta
        TABLE minf, mexp DEPEND dt FROM -100 TO 100 WITH 2000

                
        alpha = am*exp(bm*cm*(v-dm))
        beta = am*exp(-bm*(1-cm)*(v-dm))
        tau = 1/(alpha + beta)
        minf = alpha*tau
        if (tau<taum_min) { tau = taum_min }
        mexp = 1 - exp(-dt/tau)

}
 
 
UNITSON
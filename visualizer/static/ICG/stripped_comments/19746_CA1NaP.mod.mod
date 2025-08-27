UNITS {
        (mA) = (milliamp)
        (mV) = (millivolt)
}
 
NEURON {
        SUFFIX CA1NaP
        USEION na READ ena WRITE ina
        RANGE gnaP
        GLOBAL mPinf, mPexp
}
 
INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}
 
PARAMETER {
        v (mV)
        dt (ms)
        gnaP = 0.00017 (mho/cm2)
        
        vhalf = -49 (mV)
}
 
STATE {
        mP
}
 
ASSIGNED {
        ena (mV)
        ina (mA/cm2)
        mPinf mPexp
}
 
BREAKPOINT {
        SOLVE states
        ina = gnaP*mP*(v - ena)
}
 
UNITSOFF
 
INITIAL {
	rates(v)
	mP = mPinf
}

PROCEDURE states() {  
        rates(v)      
        mP = mP + mPexp*(mPinf-mP)
}
 
PROCEDURE rates(v) {  
                      
        LOCAL alpha, beta
        
                
        alpha = -1.74*(v-11)/(exp(-(v-11)/12.94)-1)
        beta = 0.06*(v-5.9)/(exp((v-5.9)/4.47)-1)
        mPinf = 1/(1+exp(-(v-vhalf)/5))
        mPexp = 1 - exp(-dt*(alpha+beta))
}
 
 
UNITSON
UNITS {
        (mA) = (milliamp)
        (mV) = (millivolt)
}
 
NEURON {
        SUFFIX KDcvode
        USEION k READ ek  WRITE ik
        RANGE  gkbar, ik, gk
} 
 
INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}
 
PARAMETER {
        v (mV)
        celsius = 37 (degC)
        mon = 1
	  hon = 1
        gkbar = .0045 (mho/cm2)
        
}
 
STATE {
        m h
}
 
ASSIGNED {
        ek (mV)
        ik (mA/cm2)
        gk minf hinf tau q10 alpha beta sum 
}
 
BREAKPOINT {
        SOLVE state METHOD cnexp
        gk = gkbar * m*h
	  ik = gk* (v-ek)
}
 
UNITSOFF
 
INITIAL {
	
	m = minf
	h = hinf
}


DERIVATIVE state {  
                      
        
        q10 = 3^((celsius - 37)/10)
        
                

        alpha = 8.5/(1+exp((v+17)/(-12.5)))
        beta =  35/(1+exp((v+99)/14.5))
        sum = alpha + beta
        minf = alpha/sum
	  tau= 10/(q10 *sum)
        m' = mon * (minf-m)/tau      

                

        alpha = 0.0015/(1+exp((v+89)/8))
        beta = 0.0055/(1+exp((v+83)/(-8)))
        sum = alpha + beta
        hinf = alpha/sum
	  tau = 1/(q10 * sum *1.6)
     	  h' = hon * (hinf-h)/tau
}

 
UNITSON
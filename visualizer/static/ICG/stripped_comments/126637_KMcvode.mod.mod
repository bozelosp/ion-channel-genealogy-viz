UNITS {
        (mA) = (milliamp)
        (mV) = (millivolt)
}
 
NEURON {
        SUFFIX KMcvode
	  USEION k READ ek WRITE ik
        RANGE  gkbar, gk, minf
} 
 
INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}
 
PARAMETER {
        v (mV)
        celsius = 37 (degC)
        mon = 1
        gkbar	= .00004 (mho/cm2)
        

}
 
STATE {
        m 
}
 
ASSIGNED {
        ek (mV)
        ik (mA/cm2)
        gk minf tau q10 alpha beta sum 
}
 
BREAKPOINT {
        SOLVE state METHOD cnexp 
        gk = gkbar *m
    	  ik = gk* (v-ek)
}
 
UNITSOFF
 
INITIAL {
	
	m = minf
}


DERIVATIVE state {  
                      

        q10 = 3^((celsius - 37)/10)
        
                

        sum = 3.3*(exp((v+35)/40)+exp(-(v+35)/20))/200
        minf = 1.0 / (1+exp(-(v+35)/10))
	  tau= 1/(q10 * sum)
        m' = mon * (minf-m)/tau      
              
}

 
UNITSON
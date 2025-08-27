UNITS {
        (mA) = (milliamp)
        (mV) = (millivolt)
}
 
NEURON {
        SUFFIX Khcvode
  	  USEION k WRITE ik
        RANGE  gkbar, gk, ik
} 
 
INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}
 
PARAMETER {
        v (mV)
        celsius = 37 (degC)
        mon = 1        
	  man = 1
	  nan = 1
        gkbar = .0003 (mho/cm2)
        ekcvode	= -30 (mV)

}
 
STATE {
        m 
}
 
ASSIGNED {
        ik (mA/cm2)
        gk minf mtau ntau q10
}
 
BREAKPOINT {
        SOLVE state METHOD cnexp
        gk = gkbar *m
	ik = gk* (v-ekcvode)
}
 
UNITSOFF
 
INITIAL {
	
	m = minf
}


 
DERIVATIVE state {  
                      

        q10 = 3^((celsius - 37)/10)
        	
                

        minf = 1/(1+exp((v+78)/7))
	  mtau= 38/q10
        ntau = 319/q10
        m' = mon * ( (man * 0.8 * (minf-m)/mtau) + (nan * 0.2 * (minf-m)/ntau) )     
}

 
UNITSON
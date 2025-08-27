NEURON {
	SUFFIX hcb
	NONSPECIFIC_CURRENT i
	
    RANGE  gbar,vhalf, K, taun, ninf, g, ihi 
	GLOBAL eh

}

UNITS {
	(um) = (micrometer)
	(mA) = (milliamp)
	(uA) = (microamp)
	(mV) = (millivolt)
	(pmho) = (picomho)
	(mmho) = (millimho)
}



PARAMETER {              
        ena    = 55    (mV)
        eh     = -10   (mV)
	K      = 10.0   (mV)	
	gbar   = 1     (mho/cm2)  
	vhalf  = -90   (mV)       
}	


STATE {                
	n
}

ASSIGNED {            
        v 
	i (mA/cm2)
	ninf
	taun (ms)
	g
}

        


INITIAL {               
	rates()	
	n = ninf
	g = gbar*n
	i = g*(v-eh)
}


BREAKPOINT {
	SOLVE states METHOD cnexp
	g = gbar*n 
	i = g*(v-eh)  
}

DERIVATIVE states {
	rates()
        n' = (ninf - n)/taun
}

PROCEDURE rates() {  
 
 	if (v > -10) {
	   taun = 1
	} else {
           taun = 2*(1/(exp((v+145)/-17.5)+exp((v+16.8)/16.5)) + 10) 

	}  
         ninf = 1 - (1 / (1 + exp((vhalf - v)/K)))                  
}
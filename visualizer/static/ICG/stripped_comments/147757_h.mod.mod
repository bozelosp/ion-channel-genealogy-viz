NEURON {
	SUFFIX h
        RANGE  gbar,vhalf, K, taun, ninf, g  
	
	NONSPECIFIC_CURRENT i
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
        dt             (ms)
	v              (mV)
        
        
	K      = 8.5   (mV)
	gbar   = 1.0     (mho/cm2)  
	
      vhalf  = -81   (mV)       
     

}	


STATE {                
	n
}

ASSIGNED {             
	i (mA/cm2)
	ninf
	taun (ms)
	g
        eh (mV)
}

        


INITIAL {               
	states()	
	n = ninf
	g = gbar*n
	i = g*(v-eh)
}


BREAKPOINT {
	SOLVE h METHOD derivimplicit
	g = gbar*n
	i = g*(v-eh)  
}

DERIVATIVE h {
	states()
        n' = (ninf - n)/taun
}

PROCEDURE states() {  
 
 	if (v > -30) {
	   taun = 1
	} else {
           
           taun = 5*(1/(exp((v+145)/-17.5)+exp((v+16.8)/16.5)) + 5) 


	}  
         ninf = 1 - (1 / (1 + exp((vhalf - v)/K)))                  
}
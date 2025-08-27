NEURON {
	SUFFIX h
        RANGE  gbar,vhalf, K, taun, ninf, g ,ina
	USEION na READ ena WRITE ina      

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
        ena    = 50    (mV)
        eh     = -10   (mV)
	K      = 8.5   (mV)
	gbar   = 0     (mho/cm2)  
	
      vhalf  = -81   (mV)       
     

}	


STATE {                
	n
}

ASSIGNED {             
	ina (mA/cm2)
	ninf
	taun (ms)
	g
}

        


INITIAL {               
	states()	
	n = ninf
	g = gbar*n
	ina = g*(v-eh)
}


BREAKPOINT {
	SOLVE h METHOD derivimplicit
	g = gbar*n
	ina = g*(v-eh)  
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
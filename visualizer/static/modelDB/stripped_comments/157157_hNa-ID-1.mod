NEURON {
	SUFFIX hNa
        RANGE  gbar,vhalf, K, taun, ninf, g  
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

INDEPENDENT {t FROM 0 TO 1 WITH 100 (ms)}

PARAMETER {              
        dt             (ms)
	v              (mV)
        ena    = 50    (mV)
        eh     = -10   (mV)
	K      = 8.5   (mV)

	gbar   = 0     (mho/cm2)  
	vhalf  = -90   (mV)       
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

	ina = g*(v-eh)*0.001            
}


BREAKPOINT {
	SOLVE h METHOD derivimplicit
	g = gbar*n

	ina = g*(v-eh)*0.001            
}

DERIVATIVE h {
	states()
        n' = (ninf - n)/taun
}

PROCEDURE states() {  
 
 	if (v > -30) {
	   taun = 1
	} else {
           taun = 2*(1/(exp((v+145)/-17.5)+exp((v+16.8)/16.5)) + 5) 

	}  
         ninf = 1 - (1 / (1 + exp((vhalf - v)/K)))                  
}
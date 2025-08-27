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

INDEPENDENT {t FROM 0 TO 1 WITH 100 (ms)}

PARAMETER {                    
        dt             (ms)
	v              (mV)
        

	K      = 8.5   (mV)
	gbar   = 0.1   (mmho/cm2) 
	vhalf  =-90    (mV)
}	


STATE {                 
	n
}

ASSIGNED {             
	eh (mV)
        i (mA/cm2)
	ninf
	taun (ms)
	g
}

        


INITIAL {      
	states()	
	n = ninf
	g = gbar*n
	i = g*(v-eh)*(0.001)
}


BREAKPOINT {
	SOLVE h METHOD derivimplicit
	g = gbar*n
	i = g*(v-eh)*(0.001)
}

DERIVATIVE h {
	states()
        n' =(ninf - n)/taun
}

PROCEDURE states() {  
 
 	if (v > -30) {
	   taun = 1
	} else {

            taun = -0.5(ms/mV)*v - 10(ms)
	}
 
        ninf = 1 - (1 / (1 + exp((vhalf - v)/K)))

}
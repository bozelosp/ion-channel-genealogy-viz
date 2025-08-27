NEURON {
	SUFFIX h
	USEION na READ ena WRITE ina      
        RANGE  gbar,vhalf, K, taun, ninf,ina

}

UNITS {
	(um) = (micrometer)
	(mA) = (milliamp)
	(uA) = (microamp)
	(mV) = (millivolt)
}

PARAMETER {              
        eh     = -10   (mV)
	K      = 8.5   (mV)
	gbar   = 0     (mho/cm2)  
	vhalf  = -90   (mV)       
}	

ASSIGNED {             
	v              (mV)
        ena            (mV)
	ina            (mA/cm2)
	ninf
	taun           (ms)
}

STATE {                
	n
} 


INITIAL {               
	states()	
	n = ninf
}


BREAKPOINT {
	SOLVE states METHOD cnexp
	ina = gbar*n*(v-eh)            
}

DERIVATIVE states {
	rates(v)
        n' = (ninf - n)/taun
}

PROCEDURE rates(v (mV)) {  
 
 	if (v > -30) {
	   taun = 1(ms)
	} else {
           taun = 2(ms)*(1/(exp((v+145)/-17.5(mV))+exp((v+16.8)/16.5(mV))) + 5) 

	}  
         ninf = 1 - (1 / (1 + exp((vhalf - v)/K)))                  
}
NEURON {
	  SUFFIX nap
	  USEION na READ ena WRITE ina
    RANGE  gnabar,vhalf, K, g, gmax
}

UNITS {
	  (mA) = (milliamp)
	  (mV) = (millivolt)
}

PARAMETER { 
	  v               (mV)
    ena = 50        (mV) 
	  K = 1           (mV)              
    
	  gnabar = 0      (mho/cm2)          
	  vhalf  = -51.90 (mV)                
}	

STATE { n }

ASSIGNED {
	  ina  (mA/cm2)
    g    (mho/cm2)
    gmax (mho/cm2)
}

INITIAL {
    
    gmax = 0
}

BREAKPOINT {
    states(v)
    g = gnabar*n*n*n
	  ina = g*(v-ena)
    if (g > gmax) {
        gmax = g
    }
}

PROCEDURE states(v (mV)) {     
    
    TABLE n DEPEND vhalf, K FROM -150 TO 150 WITH 300     
    n = 1 / (1 + exp((vhalf - v)/K)) 
}
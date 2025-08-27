NEURON {
	SUFFIX nap
	USEION na READ ena WRITE ina
        RANGE  gnabar,vhalf, K

}

UNITS {
	(mA) = (milliamp)
	(mV) = (millivolt)

}

INDEPENDENT {t FROM 0 TO 1 WITH 100 (ms)}

PARAMETER {                         
        dt (ms) 
	v (mV)
        ena = 50           (mV)     
	K = 4.5            (1)      

	gnabar = 1.0                  
	vhalf  = -50.4    (mV)      
}	

STATE { n }

ASSIGNED {
	ina (mA/cm2)
}

INITIAL {

}

BREAKPOINT {
	SOLVE states
	ina = gnabar*n*n*n*(v-ena)
}

PROCEDURE states() {     



        n = 1 / (1 + (exp(vhalf - v)/K)) 
        VERBATIM
        return 0;
        ENDVERBATIM
}
NEURON {
	SUFFIX nap
	USEION na READ ena WRITE ina
        RANGE  gnabar,vhalf, K, ina

}

UNITS {
	(mA) = (milliamp)
	(mV) = (millivolt)

}

PARAMETER {            
	K = 4.5            (1)      

	gnabar = 0        (mho/cm2)
	vhalf  = -50.4    (mV)      
      
}	

ASSIGNED {
	v             (mV)
        ena           (mV)    
	ina           (mA/cm2)
        n_inf
        tau            (ms)
}

STATE { n }

BREAKPOINT {
	SOLVE states METHOD cnexp
	ina = gnabar*n*n*n*(v-ena)
}

INITIAL {
	rates(v)
	n = n_inf
}


DERIVATIVE states {
        rates(v)
        n' = (n_inf-n)/tau
}

PROCEDURE rates(v(mV)) {
	n_inf = 1 / (1 + (exp(vhalf - v)/K))
	tau =1
}
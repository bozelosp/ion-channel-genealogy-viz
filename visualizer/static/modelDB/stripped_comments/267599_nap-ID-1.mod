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
	K = 2            (1)      

	gnabar = 0        (mho/cm2)
	vhalf  = -60.4    (mV)      
      
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
    LOCAL Arg
    tau=10
    Arg=(vhalf-v)/K

	if (Arg<-50) {n_inf=1}
    else if (Arg>50) {n_inf=0}
    else {n_inf=1/(1+exp(Arg))}
}
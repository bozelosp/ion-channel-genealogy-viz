NEURON {
	SUFFIX IK1
	USEION k READ ek WRITE ik
	RANGE gK1, ik
	GLOBAL dummy 

}

UNITS {
	(mA) = (milliamp)
	(mV) = (millivolt)
        (mM) = (milli/liter)
	 

}

PARAMETER {
	 gK1=0.1567e-3 (S/cm2) <0,1e9>
	
}


ASSIGNED {
	v (mV)
    	ik (mA/cm2)
	ek (mV)
	dummy
	      
}

LOCAL k
BREAKPOINT {
	
	ik = gK1/(1 + exp(0.07*(v + 80)))*(v - ek)
}
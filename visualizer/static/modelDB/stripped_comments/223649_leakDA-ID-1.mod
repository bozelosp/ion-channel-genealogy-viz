NEURON {
	SUFFIX leakDA
	NONSPECIFIC_CURRENT il
	RANGE il, el, glbar
}

UNITS {
	(mA) = (milliamp)
	(mV) = (millivolt)
}

PARAMETER {
	DA_start = 100		             
	DA_stop = 600	
	DA_t1 = 0.8 
	
	glbar = 2.857142857142857e-05  
	el = -75 (mV)
}

ASSIGNED {
	v (mV)
	il (mA/cm2)
}

BREAKPOINT { 
	il = glbar*(v - el)*DA1(t)
}
FUNCTION DA1(t) {
	    if (t >= DA_start && t <= DA_stop){DA1 = DA_t1} 									
		else  {DA1 = 1}
	}
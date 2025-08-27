NEURON { 
	SUFFIX GRC_LKG3 
	NONSPECIFIC_CURRENT il
	RANGE el, gl,i
} 
 
UNITS { 
	(mA) = (milliamp) 
	(mV) = (millivolt) 
} 
 
PARAMETER { 
	v (mV) 
	gl = 5e-7 (mho/cm2)
	celsius = 30 (degC)
	el =  -70 (mV) 
	
} 

ASSIGNED { 
	il (mA/cm2) 
	i (mA/cm2) 
}
  
BREAKPOINT { 
	il = gl*(v - el) 
	i = il
}
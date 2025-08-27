NEURON { 
	SUFFIX GrC_Lkg1 
	NONSPECIFIC_CURRENT il
	RANGE el, gl,i
} 
 
UNITS { 
	(mA) = (milliamp) 
	(mV) = (millivolt) 
} 
 
PARAMETER { 
		gl = 5.68e-5 (mho/cm2)
	      el = -58 (mV)
} 

ASSIGNED {
      v (mV) 
      il (mA/cm2) 
	i (mA/cm2) 
}
  
BREAKPOINT { 
	il = gl*(v - el) 
	i = il
}
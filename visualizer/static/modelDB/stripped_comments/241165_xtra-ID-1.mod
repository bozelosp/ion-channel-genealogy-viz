NEURON {
	SUFFIX xtra
	RANGE es 
	RANGE x, y, z, type, order
	GLOBAL stim 
	POINTER ex 
}

PARAMETER {	
	es = 0 (mV)
	x = 0 (1) 
	y = 0 (1)
	z = 0 (1)		
	type = 0 (1) 
	order = 0 (1) 
}

ASSIGNED {
	v (millivolts)
	ex (millivolts)
	stim (unitless) 		
	area (micron2)
}

INITIAL {
	ex = stim*es	
}


BEFORE BREAKPOINT { 
  ex = stim*es
}
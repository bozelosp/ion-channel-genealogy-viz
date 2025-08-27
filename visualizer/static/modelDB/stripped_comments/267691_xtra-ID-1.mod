NEURON {
	SUFFIX xtra
	RANGE es1 
	RANGE es2 
	RANGE ampratio	
	RANGE x, y, z, type, order
	GLOBAL stim 
	POINTER ex 
}

PARAMETER {	
	es1 = 0 (mV)
    es2 = 0 (mV)
	x = 0 (1) 
	y = 0 (1)
	z = 0 (1)		
	type = 0 (1) 
	order = 0 (1) 
    ampratio = 0	
}

ASSIGNED {
	v (millivolts)
	ex (millivolts)
	stim (unitless) 		
	area (micron2)
}

INITIAL {
	ex = (stim*es1)+(stim*ampratio*es2)	
}


BEFORE BREAKPOINT { 
  ex = (stim*es1)+(stim*ampratio*es2)
}
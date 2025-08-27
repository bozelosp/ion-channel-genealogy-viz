INCLUDE "custom_code/inc_files/144520_Unit.inc"
INCLUDE "custom_code/inc_files/144520_Volume.inc"
NEURON {
	SUFFIX ibg 
	USEION na READ ena WRITE ina
	
	RANGE gCa, gNa, ina, ica
}
PARAMETER {
	gNa = 0.18	(uS)		<0,1e9>
	
}

ASSIGNED {
	celsius (degC)
	v (mV) 
	ina (mA/cm2)  
	
	ena (mV)
	
	dummy
}

BREAKPOINT {
	ina =  (1e-06)*gNa/S*(v - ena)
	
}
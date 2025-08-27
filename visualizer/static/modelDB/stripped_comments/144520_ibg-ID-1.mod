INCLUDE "Unit.inc"
INCLUDE "Volume.inc"
NEURON { SUFFIX ibg 
	USEION na READ ena WRITE ina
	USEION ca READ eca WRITE ica
	RANGE gCa, gNa, ina, ica
}
PARAMETER {
	gNa = 0.18	(uS)		<0,1e9>
	gCa = 0.02	(uS)		<0,1e9>
}

ASSIGNED {
	celsius (degC)
	v (mV) 
	ina (mA/cm2)  
	ica (mA/cm2)  
	ena (mV)
	eca (mV)
	dummy
}

BREAKPOINT {
	ina =  (1e-06)*gNa/S*(v - ena)
	ica =  (1e-06)*gCa/S*(v - eca)
}
NEURON {
	SUFFIX hdend
	USEION na READ ena WRITE ina
	USEION k  READ ek  WRITE ik
	RANGE gbar
	RANGE ninf, ntau

	GLOBAL vhalf_n, vsteep_n, exp_n 
	GLOBAL tskew_n, tscale_n, toffset_n 

}

INCLUDE "noinact_nak_currs.inc"

INCLUDE "noinact_gate_states.inc"

INCLUDE "var_funcs.inc"
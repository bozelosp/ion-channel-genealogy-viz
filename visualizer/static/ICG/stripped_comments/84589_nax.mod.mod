NEURON {
	SUFFIX nax
	USEION na READ ena WRITE ina
	




	RANGE gbar
	RANGE minf, mtau, hinf, htau

	GLOBAL vhalf_m, vsteep_m, exp_m 
	GLOBAL tskew_m, tscale_m, toffset_m 
	
	GLOBAL vhalf_h, vsteep_h, exp_h
	GLOBAL tskew_h, tscale_h, toffset_h 
}


INCLUDE "custom_code/inc_files/84589_inact_na_currs.inc"

INCLUDE "custom_code/inc_files/84589_nax_inact_gate_states.inc"

INCLUDE "custom_code/inc_files/84589_var_funcs.inc"
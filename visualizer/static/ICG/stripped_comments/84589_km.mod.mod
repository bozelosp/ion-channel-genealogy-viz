NEURON {
	SUFFIX km
	USEION k READ ek WRITE ik
	RANGE gbar
	RANGE ninf, ntau

	GLOBAL vhalf_n, vsteep_n, exp_n 
	GLOBAL tskew_n, tscale_n, toffset_n 
}


INCLUDE "custom_code/inc_files/84589_noinact_k_currs.inc"

INCLUDE "custom_code/inc_files/84589_km_noinact_gate_states.inc"

INCLUDE "custom_code/inc_files/84589_var_funcs.inc"
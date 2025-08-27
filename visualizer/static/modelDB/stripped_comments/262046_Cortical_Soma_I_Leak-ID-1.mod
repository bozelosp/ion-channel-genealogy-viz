UNITS {
 (mV) = (millivolt)
 (mA) = (milliamp)
 (S) = (siemens)
}

NEURON {
	SUFFIX cortical_soma_i_leak
	NONSPECIFIC_CURRENT i_l			
	RANGE i_l, g_l, e_l				
}

PARAMETER {
	g_l = 0.1e-3 (S/cm2)
	e_l = -70 (mV)
}

ASSIGNED {
	v (mV)	
	i_l (mA/cm2)
}

BREAKPOINT {
	i_l = g_l*(v - e_l)
}
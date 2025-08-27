NEURON {
	SUFFIX Golgi_lkg
	NONSPECIFIC_CURRENT i
	RANGE Q10_diff,gbar_Q10, fix_celsius
	RANGE el, gbar, g	, ic
}

UNITS {
	(mA) = (milliamp)
	(mV) = (millivolt)
}

PARAMETER {
	v (mV)
	Q10_diff	= 1.5
	gbar = 21e-6 (mho/cm2)
    fix_celsius = 37 (degC)
	el = -55 (mV)
}

ASSIGNED {
	i (mA/cm2)
	gbar_Q10 (mho/cm2)
	ic
	g
}

INITIAL{
	gbar_Q10 = gbar*(Q10_diff^((fix_celsius-23)/10))
	g = gbar_Q10
}

BREAKPOINT {
	i = gbar_Q10 * (v - el )
	ic = i
}
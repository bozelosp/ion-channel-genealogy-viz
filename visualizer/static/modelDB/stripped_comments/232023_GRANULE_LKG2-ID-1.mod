NEURON {
	SUFFIX GRANULE_LKG2
	NONSPECIFIC_CURRENT il
	RANGE Q10_diff,Q10_channel,gbar_Q10, fix_celsius
	RANGE egaba, g , ic, gbar
}

UNITS {
	(mA) = (milliamp)
	(mV) = (millivolt)
}

PARAMETER {
	v (mV)
	gbar = 6e-5 (mho/cm2) 
	egaba = -65 (mV)
	Q10_diff	= 1.5
    fix_celsius = 37 (degC)
}

ASSIGNED {
	il (mA/cm2)
	ic (mA/cm2)
	g (mho/cm2)
	gbar_Q10 (mho/cm2)
}

INITIAL {
  gbar_Q10 = gbar*(Q10_diff^((fix_celsius-30)/10))
  g = gbar_Q10
}
BREAKPOINT {
    il = g*(v - egaba)
    ic =il
}
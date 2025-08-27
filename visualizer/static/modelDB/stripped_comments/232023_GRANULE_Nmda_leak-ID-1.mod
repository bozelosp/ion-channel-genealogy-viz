NEURON {
	SUFFIX GRANULE_Nmda_leak
	NONSPECIFIC_CURRENT i
	RANGE Q10_diff
	RANGE g , ic, Erev, gmax, gmax_factor, fix_celsius
	RANGE MgBlock,v0_block,k_block
}

UNITS {
    (nA) = (nanoamp)
    (mV) = (millivolt)
    (umho) = (micromho)
    (mM) = (milli/liter)
    (pS) = (picosiemens)
    (nS) = (nanosiemens)
    (um) = (micrometer)
    PI	= (pi)		(1)
}

PARAMETER {
    gmax_factor = 1
    gmax	= 50e-12	(mho)
    surf        = 299.26e-8 (cm2)
    Erev	= -3.7  (mV)	
    Q10_diff	= 1.5
    v0_block = -20 (mV)
    k_block  = 13 (mV)
    O = 10e-3 (s) 
    Or = 52 
    v		(mV)		
    fix_celsius = 37 (degC)
}

ASSIGNED {
    i 		(mA/cm2)		
    ic 		(mA/cm2)		
    g 		(mho/cm2)		
    MgBlock
    gbar_Q10 (1)
}

INITIAL {
	rates(v)
	gbar_Q10 = Q10_diff^((fix_celsius-30)/10)
}

BREAKPOINT {
	rates(v)
	g = gmax / surf * gbar_Q10 * O * Or * gmax_factor  * MgBlock
	i = g * (v - Erev)
	ic = i
}

PROCEDURE rates(v(mV)) {
	
	TABLE MgBlock DEPEND v0_block,k_block FROM -120 TO 30 WITH 150
	MgBlock = 1 / ( 1 + exp ( - ( v - v0_block ) / k_block ) )
}
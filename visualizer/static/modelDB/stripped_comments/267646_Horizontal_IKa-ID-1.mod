INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}

NEURON {
	SUFFIX Hor_IKa
	USEION k READ ek WRITE ik
	RANGE gbar
}

UNITS {
	(mA) = (milliamp)
	(mV) = (millivolt)
}

PARAMETER {
    gbar = 4.2453e-05 (mho/cm2) 
}

ASSIGNED {
    v    (mV)
    ek   (mV)
	ik   (mA/cm2)

    m    (mV)
}

BREAKPOINT {
UNITSOFF
    m = 1 / ( 1 + exp( (v+60)/12 ) )
UNITSON
    ik = gbar * m*m*m*m*m * (v - ek)
}
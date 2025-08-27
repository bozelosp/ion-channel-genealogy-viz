INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}

NEURON {
	SUFFIX kl
        USEION k READ ek WRITE ik VALENCE 1
	RANGE gmax, i
        GLOBAL erev
}

UNITS {
	(nA) = (nanoamp)
	(mV) = (millivolt)
	(umho) = (micromho)
}

PARAMETER {
  ek  (mV)
  gmax	= 4e-4	(S/cm2)		
  erev	= -100	(mV)		
}

ASSIGNED {
	v		(mV)		
	i 		(mA/cm^2)
	ik 		(mA/cm^2)
}

INITIAL {
}

BREAKPOINT {
	i = gmax * (v - ek)
        ik=i
}
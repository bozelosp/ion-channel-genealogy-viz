INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}

NEURON {
	POINT_PROCESS noise
	NONSPECIFIC_CURRENT i
	RANGE imax

}

ASSIGNED {
	rn
}

UNITS {
	(nA) = (nanoamp)
}

PARAMETER {
	imax=1	        (umho)

}

INITIAL {

	rn = (1/(2^31))*2
}

ASSIGNED { i (nA) }

BREAKPOINT {

SOLVE dum	

}

PROCEDURE dum() {
i = (scop_random()*rn-1)*imax
}
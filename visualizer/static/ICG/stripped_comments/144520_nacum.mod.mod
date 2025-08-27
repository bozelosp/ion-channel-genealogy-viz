INCLUDE "Unit.inc"
INCLUDE "Volume.inc"
NEURON {
	SUFFIX Na_acc
	USEION na READ ina, nai  WRITE nai
	RANGE ina
}

PARAMETER {
	ina			(mA/cm2)
}

STATE {
	nai START 8	(mM)
}

LOCAL ViF
INITIAL {
	VERBATIM
		nai = _ion_nai;
	ENDVERBATIM
	ViF = (1e-3)*Vi*F/S
}

BREAKPOINT {
	SOLVE state METHOD derivimplicit
}

DERIVATIVE state {
	nai' = -ina/(ViF)
}
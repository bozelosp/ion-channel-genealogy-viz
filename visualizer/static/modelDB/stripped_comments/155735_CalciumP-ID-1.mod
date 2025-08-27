INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}

NEURON {
	SUFFIX cad
	USEION ca READ ica, cai WRITE cai
	RANGE ca, depth
	GLOBAL cainf,taur
}

UNITS {
	(molar) = (1/liter)      
	(mM)	= (millimolar)
	(um)	= (micron)
	(mA)	= (milliamp)
	(msM)   = (ms mM)
}

PARAMETER {
	depth = .2	(um)     
	taur  = 0.1     (ms)     
	cainf = 4e-5	(mM)
	cai		(mM)
}

STATE {
	ca   (mM)
}

INITIAL {
	ca = cainf
	cai = ca
}

ASSIGNED{
	ica		(mA/cm2)
	drive_channel   (mM/ms)
}

BREAKPOINT{
	SOLVE state METHOD cnexp
}

DERIVATIVE state {

	drive_channel = -(10000)*ica/(2*96494*depth)

	if(drive_channel <= 0.) {drive_channel = 0.}

	ca' = drive_channel + (cainf-ca)/taur
	cai = ca
}
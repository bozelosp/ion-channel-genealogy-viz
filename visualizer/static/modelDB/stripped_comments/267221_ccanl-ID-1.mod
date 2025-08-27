NEURON {
	SUFFIX ccanl
	USEION ca READ cai, ica, eca WRITE eca, cai
	RANGE caiinf, catau, cai, eca
}

UNITS {
	(molar) = (1/liter)
	(mM) = (millimolar)
	(um) = (micron)
	(mV) = (millivolt)
	(mA) = (milliamp)
	FARADAY = (faraday) (coulomb)
	R =(k-mole) (joule/degC)
}

PARAMETER {
	celsius = 6.3 (degC)
	depth = .2 (um)	
	catau = 9 (ms)

	
	
	caiinf = 5e-5 (mM)
	
	cao = 2 (mM)
}

ASSIGNED {
	ica (mA/cm2)
	eca (mV)
	drive_channel (mM/ms)
	cai (mM)
}

STATE {
	ca (mM)
}

BREAKPOINT {
	SOLVE states METHOD euler
	eca = KTF(celsius) * log(cao/ca)
}

DERIVATIVE states {

	drive_channel =  -(1e4)*ica/(depth*FARADAY)
	if (drive_channel <= 0.) { drive_channel = 0.  }

	ca' = drive_channel + (caiinf/3 - ca)/catau
	cai = ca
}

INITIAL {
	ca = caiinf	
	eca = KTF(celsius) * log(cao/ca)
}

FUNCTION KTF(celsius (degC)) (mV) {
	KTF = (1e3)*R*(celsius + 273.15(degC))/(2*FARADAY)
}
NEURON {
	SUFFIX nadifl
	USEION na READ ina WRITE nai
	RANGE D
}

UNITS {
	(mM) = (milli/liter)
	(um) = (micron)
	FARADAY = (faraday) (coulomb)
	PI = (pi) (1)
}

PARAMETER {
	D = .0 (um2/ms)
}

ASSIGNED {
	ina (milliamp/cm2)
	diam (um)
}

STATE {
	nai (mM)
}

BREAKPOINT {
	SOLVE conc METHOD sparse
}

KINETIC conc {
	COMPARTMENT PI*diam*diam/4 {nai}
	LONGITUDINAL_DIFFUSION D*PI*diam*diam/4 {nai}
	~ nai << (-ina/(FARADAY)*PI*diam*(1e4))
}
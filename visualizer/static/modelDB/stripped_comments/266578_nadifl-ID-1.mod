NEURON {
	SUFFIX nadifl
	USEION na READ ina WRITE nai, ena

	RANGE D, Nai, Total, neo, tau
}

UNITS {
	(mM) = (milli/liter)
	(um) = (micron)
	FARADAY = (faraday) (coulomb)
	PI = (pi) (1)
}

PARAMETER {
	D = .6 (um2/ms)
         k1buf = 0.01          
        k2buf =  0          



       tau = 5000 (ms)
}

ASSIGNED {
	ina (milliamp/cm2)
	diam (um)
        iNa (milliamp/cm2)
        Total (mM)
        neo (milliamp/cm2)
        ena (mV)
}

STATE {
	nai (mM)
        Nai (mM)
}

INITIAL {
lates()

nai = 10
Nai = 10
Total = 20
neo = ina	
}

BREAKPOINT {
	SOLVE conc METHOD sparse
}

KINETIC conc {
        lates()
	COMPARTMENT PI*diam*diam/4 {nai}
	LONGITUDINAL_DIFFUSION D {nai}
~ nai << (-neo/(FARADAY)*PI*diam*(1e4))







}


PROCEDURE lates() {
LAG ina BY tau
  neo = lag_ina_tau
if (ena < 70) {ena = 70}
if (nai < 10) {nai = 10}
}
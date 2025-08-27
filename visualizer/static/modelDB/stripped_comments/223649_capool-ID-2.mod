NEURON {
    SUFFIX capool
	USEION ca READ ica
	USEION cas READ casi WRITE casi VALENCE 2 
	RANGE fcas, taucas, cainf
}

UNITS {
        (mM) = (milli/liter)
        (mA) = (milliamp)
	(mV) = (millivolt)
	FARADAY = 96485.309 (coul)
}

PARAMETER {
	pi = 3.14159265
	taucas= 1000 (ms) 	
    cainf= 50e-6   (mM)  	
	fcas = 0.024
    w = 1 (micrometer)     	
	z = 2			
}

ASSIGNED {
	v (mV)
	ica (mA/cm2)
    A       (mM-cm2/ms/mA)
}

STATE { casi(mM) }

BREAKPOINT { 
	SOLVE states METHOD cnexp
}

INITIAL {
	A = 1/(z*FARADAY*w)*(1e4)
	casi = cainf
}

DERIVATIVE states {
	casi'= -fcas*A*ica + (cainf - casi)/taucas
}
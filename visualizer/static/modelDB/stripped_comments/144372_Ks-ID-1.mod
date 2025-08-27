NEURON {
        SUFFIX Ks
        USEION k READ ek WRITE ik
        RANGE g, ik
	GLOBAL tau
}

UNITS {
	(S)  = (siemens)
        (mA) = (milliamp)
        (mV) = (millivolt)
}

PARAMETER {
        g 	(S/cm2)
	tau= 75 (ms)
}

ASSIGNED {
        v 	(mV)
	ek 	(mV)
        ik 	(mA/cm2)
	mtau 	(ms)
	minf
}

STATE { m }

BREAKPOINT { 
        SOLVE states METHOD cnexp
	ik= g* m* (v- ek) 
}

DERIVATIVE states {
	rates()
	m'= (minf- m)/ mtau
}

INITIAL {
	rates()
	m= minf 
}

PROCEDURE rates() { UNITSOFF
	minf= 1/ (1+ exp(-(v+ 39)/ 5))
	mtau= tau
} UNITSON
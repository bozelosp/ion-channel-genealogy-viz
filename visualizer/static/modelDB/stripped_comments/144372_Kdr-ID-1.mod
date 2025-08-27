NEURON {
        SUFFIX Kdr
        USEION k READ ek WRITE ik
        RANGE g, ik
}

UNITS {
	(S)  = (siemens)
        (mA) = (milliamp)
        (mV) = (millivolt)
}

PARAMETER {
        g	(S/cm2)
}

ASSIGNED {
        v	(mV)
	ek	(mV)
        ik	(mA/cm2)
	ntau	(ms)
	ninf
}

STATE { n }

BREAKPOINT { 
        SOLVE states METHOD cnexp
	ik= g* n^4* (v- ek) 
}

DERIVATIVE states {
	rates()
	n'= (ninf- n)/ ntau
}

INITIAL {
	rates()
	n= ninf 
}

PROCEDURE rates() { UNITSOFF
	ninf= 1/ (1+ exp(-(v+ 30)/ 10))
	ntau= 0.37+ 1.85/ (1+ exp((v+ 27)/ 15))
} UNITSON
NEURON {
        SUFFIX H
	NONSPECIFIC_CURRENT i
        RANGE g, i, hinf, h1tau, h2tau
}

UNITS {
	(S)  = (siemens)
        (mA) = (milliamp)
        (mV) = (millivolt)
}

PARAMETER {
	e= -43	(mV)
        g	(S/cm2)
}

ASSIGNED {
        v	(mV)
        i	(mA/cm2)
	h1tau	(ms)
	h2tau	(ms)
	hinf
}

STATE { h1 h2 }

BREAKPOINT { 
        SOLVE states METHOD cnexp
        i= g* (0.8* h1+ 0.2* h2)* (v- e) 
}

DERIVATIVE states {
        rates()
        h1'= (hinf- h1)/ h1tau
        h2'= (hinf- h2)/ h2tau
}

INITIAL {
	rates()
	h1= hinf
	h2= hinf
}

PROCEDURE rates() { UNITSOFF
	hinf= 1/ (1+ exp((v+ 82)/ 7))
	h1tau= 40
	h2tau= 300
} UNITSON
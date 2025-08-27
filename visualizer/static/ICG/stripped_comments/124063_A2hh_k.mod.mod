NEURON {
	SUFFIX HHk
	USEION k READ ek WRITE ik
	RANGE gkbar, ik
	GLOBAL inf, tau
}

UNITS {
	(mA) = (milliamp)
	(mV) = (millivolt)
}

PARAMETER {
	gkbar=.00026 (mho/cm2) <0,1e9>
	ek (mV) 
}
STATE {
	n
}
ASSIGNED {
	v (mV)
	celsius (degC) 
	ik (mA/cm2)
	inf
	tau (ms)
}

INITIAL {
	rate(v*1(/mV))
	n = inf
}

BREAKPOINT {
	SOLVE states METHOD cnexp
	ik = gkbar*n*n*n*n*(v - ek)
}

DERIVATIVE states {	
	rate(v*1(/mV))
	n' = (inf - n)/tau
}

UNITSOFF
FUNCTION alp(v(mV)) { LOCAL q10
	v = -v - 65
	q10 = 3 
	alp = q10 * .01*expM1(v + 10, 10)
}

FUNCTION bet(v(mV)) { LOCAL q10
	v = -v - 65
	q10 = 3^((celsius - 6.3)/10)
	bet = q10 * .125*exp(v/80)
}

FUNCTION expM1(x,y) {
        if (fabs(x/y) < 1e-6) {
                expM1 = y*(1 - x/y/2)
        }else{
                expM1 = x/(exp(x/y) - 1)
        }
}


PROCEDURE rate(v) {LOCAL a, b 
	TABLE inf, tau DEPEND celsius FROM -100 TO 100 WITH 200
		a = alp(v)/5  b=bet(v)/5
		tau = 1/(a + b)
		inf = a/(a + b)
}
UNITSON
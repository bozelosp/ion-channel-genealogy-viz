NEURON {
	SUFFIX sk
	USEION k READ ek WRITE ik
        USEION ca READ cai
        RANGE  gbar,gkahp,ik, inf,tau,g
        GLOBAL Cq10
}

UNITS {
	(mA) = (milliamp)
	(mV) = (millivolt)
	(molar) = (1/liter)
	(mM) = (millimolar)
	(pS) = (picosiemens)
	(um) = (micron)
}

PARAMETER {
	gbar = 1.0	(pS/um2)
        n = 4
        cai = 50.e-6	(mM)
        a0 = 1.3e4	(1/ms-mM-mM-mM-mM)	
        b0 = 0.06	(1/ms)			
	    celsius = 37(degC)
	Cq10 = 3
}

STATE {	w }

ASSIGNED {
	ik	(mA/cm2)
        g	(pS/um2)
        inf
        tau	(ms)
	a	(1/ms)
        v	(mV)
        ek	(mV)
}

BREAKPOINT {
	SOLVE state METHOD cnexp
	g = gbar*w
	ik = (1e-4)* g*(v-ek)
}

INITIAL {
	rate(cai)
	w=inf
}

DERIVATIVE state {
	rate(cai)
	w' = (inf - w)/tau
}

PROCEDURE rate(cai (mM)) {
	LOCAL q10
	q10 = Cq10^((celsius - 22 (degC))/10 (degC) )
	a = a0*cai^4
	tau = q10/(a + b0)
	inf = a/(a + b0)
}
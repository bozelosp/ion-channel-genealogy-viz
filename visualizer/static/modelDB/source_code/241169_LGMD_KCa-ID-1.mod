TITLE KCa channels for LGMD SFA

UNITS {
	(mV) = (millivolt)
	(mA) = (milliamp)
	(S) = (siemens)
	(molar) = (1/liter)
	(mM) = (millimolar)
}

NEURON {
	THREADSAFE
    SUFFIX KCa

    USEION ca READ cai
    USEION k READ ek WRITE ik
    RANGE gmax, g
    GLOBAL pwr, tau, kD_ca
}

PARAMETER {
    gmax = 0.01	(S/cm2)
    kD_ca =1e-4	(mM)
    tau = 0.1	(ms)
    minca =1.2e-4 (mM)
    pwr = 3		(1)
}

ASSIGNED {
    v (mV)
    cai (mM)
    ek (mV)
    
    ik (mA/cm2)
    g  (S/cm2)
}

STATE {
    n
}

INITIAL {
    n = (cai-minca)/(cai-minca+kD_ca)
}

BREAKPOINT {
	SOLVE state METHOD derivimplicit
	if (cai < minca) {
		n=0
	}
    g = n^pwr*gmax
    ik = g*(v-ek)
}

DERIVATIVE state {
	n' = ((cai-minca)/(cai-minca+kD_ca) - n)/tau
}










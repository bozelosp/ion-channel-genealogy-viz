NEURON {
	SUFFIX dr
	
	USEION k READ ek WRITE ik
	RANGE gd, gr
	RANGE td, tr
	GLOBAL ek
}


UNITS {
	(mA) = (milliamp)
	(mV) = (millivolt)
}

PARAMETER {
	gd  = 0.0025  (mho/cm2)
	gr  = 0.0028  (mho/cm2)
	V0d = -25      (mV)
    V0r = -25       (mV)
    sd = 0.25       (/mV)
    sr = 0.25       (/mV)
    tr = 1.5         (ms)
    ena = 50 (mV)
}

STATE {
	ad ar
}

ASSIGNED {
	celsius	(degC)
	
	ik      (mA/cm2)
    v       (mV)
	
    ek      (mV)
    rho     (1)
    arinf
}

INITIAL {
    rho = 1.3^((celsius - 25 (degC))/10(degC))
    ar = 1/(1+exp(-sr*(v - V0r)))
    ad = 1/(1+exp(-sd*(v - V0d)))
}

BREAKPOINT {
	SOLVE states METHOD cnexp
    rho = 1.3^((celsius - 25 (degC))/10(degC))
	ad = 1/(1+exp(-sd*(v - V0d)))
    
	ik  = rho * gr * ar * (v - ek)
}

DERIVATIVE states {
    LOCAL phi
    phi = 3^((celsius - 25 (degC))/ 10 (degC))
    arinf = 1/(1+exp(-sr*(v - V0r)))
    ar' = phi*(arinf - ar)/tr
}
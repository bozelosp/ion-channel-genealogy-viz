: Six state kinetic P/Q-type calcium channel gating scheme
: Ref: Li L, Bischofberger J, Jonas P. 2007
: Differntial gating and recruitment of P/Q-, N-, and R-type
: Ca2+ channels in hippocampal mossy fiber boutons.
: J Neurosci 27:13420-429

NEURON {
    SUFFIX mfbpqca
    USEION ca READ eca WRITE ica
    RANGE gca, gcabar
}

UNITS {
    (mA) = (milliamp)
    (mV) = (millivolt)
    (S) = (siemens)
}

PARAMETER {
    gcabar = 0.0032 (S/cm2) 
    a1o = 5.89 (/ms)    
    b1o = 14.99 (/ms)
    V1   = 62.61(mV)
    a2o = 9.21 (/ms)    
    b2o = 6.63 (/ms)
    V2   = 33.92(mV)
    a3o = 5.20 (/ms)   
    b3o = 132.80 (/ms)
    V3   = 135.08(mV)
    a4o = 1823.18 (/ms)    
    b4o = 248.58 (/ms)
    V4   = 20.86(mV)
	a5o = 247.71 (/ms)    
    b5o = 8.28 (/ms)
}

ASSIGNED {
    v    (mV)
	eca  (mV)
    gca  (S/cm2)
    ica  (mA/cm2)
	a1   (/ms)
	b1   (/ms)
	a2   (/ms)
	b2   (/ms)
	a3   (/ms)
	b3   (/ms)
	a4   (/ms)
	b4   (/ms)
	a5   (/ms)
	b5   (/ms)
}

STATE { c0 c1 c2 c3 c4 o }

BREAKPOINT {
    SOLVE kin METHOD sparse
    gca = gcabar*o
    ica = gca*(v - eca)*(1e-3)
}

INITIAL { SOLVE kin STEADYSTATE sparse }

KINETIC kin {
    rates(v)
    ~ c0 <-> c1 (a1, b1)
    ~ c1 <-> c2 (a2, b2)
    ~ c2 <-> c3 (a3, b3)
    ~ c3 <-> c4 (a4, b4)
	~ c4 <-> o  (a5, b5)
    CONSERVE c0 + c1 + c2 + c3 + c4 + o = 1
}

PROCEDURE rates(v(millivolt)) {

	a1 = a1o*exp( v/V1)
    b1 = b1o*exp(-v/V1)
    a2 = a2o*exp( v/V2)
    b2 = b2o*exp(-v/V2)
    a3 = a3o*exp( v/V3)
    b3 = b3o*exp(-v/V3) 
    a4 = a4o*exp( v/V4)
    b4 = b4o*exp(-v/V4) 
    a5 = a5o
    b5 = b5o
}







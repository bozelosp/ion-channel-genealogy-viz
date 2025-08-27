: Six state kinetic R-type calcium channel gating scheme
: Ref: Li L, Bischofberger J, Jonas P. 2007
: Differntial gating and recruitment of P/Q-, N-, and R-type
: Ca2+ channels in hippocampal mossy fiber boutons.
: J Neurosci 27:13420-429

NEURON {
    SUFFIX mfbrca
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
    a1o = 9911.36 (/ms)    
    b1o = 60.62 (/ms)
    V1   = 67.75(mV)
    a2o = 4.88 (/ms)    
    b2o = 21.91 (/ms)
    V2   = 50.94(mV)
    a3o = 4.00 (/ms)   
    b3o = 51.30 (/ms)
    V3   = 173.29(mV)
    a4o = 256.41 (/ms)    
    b4o = 116.97 (/ms)
    V4   = 16.92(mV)
	a5o = 228.83 (/ms)    
    b5o = 1.78 (/ms)
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







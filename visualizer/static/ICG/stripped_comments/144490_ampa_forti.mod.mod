NEURON {
	  POINT_PROCESS AmpaSyn
	  RANGE  e, g, i, ica
	  NONSPECIFIC_CURRENT i
    USEION ca READ cai,cao WRITE ica
    GLOBAL total, taurise, taufast, tauslow, taurise_exp, taufast_exp, tauslow_exp, afast, aslow, normfac, T_exp, fracca
}

INCLUDE "units.inc"

PARAMETER {
	  taurise_exp =.2 (ms) <1e-9,1e9>       
	  taufast_exp = 0.61 (ms) <1e-9,1e9>    
	  tauslow_exp = 2.55 (ms) <1e-9,1e9>    
    afast = 0.75 <0,1>
	  e=0	(mV)
    fracca= 0.002        
    z = 2
    celsius = 22	(degC)
    T_exp = 22    (degC)
    q10 = 3
}

ASSIGNED {
    v       (mV)
    i       (nA)
    ica	    (nA) 	
    g       (uS)
    aslow 
    total   (uS)
    cai     (mM)
    cao     (mM)
    taurise (ms)
    taufast (ms)
    tauslow (ms)
    normfac 
}

BREAKPOINT {
	  SOLVE state METHOD cnexp
    g = (B + C - A) 
    i =   g * (1-fracca) * (v - e)
    ica = g * fracca     * ghkg(v,cai,cao,z)
}

INCLUDE "triexpsyn.inc"

INCLUDE "ghk.inc"
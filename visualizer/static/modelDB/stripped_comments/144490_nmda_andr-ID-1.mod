NEURON {
    POINT_PROCESS NmdaSyn
    RANGE  e, g, i, b, ica
    NONSPECIFIC_CURRENT i
    USEION ca READ cai,cao WRITE ica
    GLOBAL total, mg, q10, taurise, taufast, tauslow, taurise_exp, taufast_exp, tauslow_exp, afast, aslow, normfac, T_exp, K0, delta, fracca
}

INCLUDE "units.inc"

PARAMETER {
    
    taurise_exp =    6.46 (ms) <1e-9,1e9>    
    taufast_exp =  252.5  (ms) <1e-9,1e9>    
    tauslow_exp = 1660    (ms) <1e-9,1e9>    
    afast = 0.61 <0,1>
    e=0	(mV)
    mg	= 1    (mM)		
    fracca= 0.13        
    z = 2
    celsius = 22	(degC)
    T_exp = 22    (degC)
    q10 = 3       
    K0 = 4.1 (mM) 
    delta = 0.8   
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
    b
}

BREAKPOINT {
    SOLVE state METHOD cnexp
    b = mgblock(v)		
    g = (B + C - A) * b
    i =   g * (1-fracca) * (v - e)
    ica = g * fracca     * ghkg(v,cai,cao,z)
}

INCLUDE "triexpsyn.inc"

INCLUDE "ghk.inc"

FUNCTION mgblock(v(mV)) {
    TABLE 
    DEPEND mg, K0, delta, celsius
    FROM -140 TO 80 WITH 1000
    
    mgblock = 1/(1+(mg/K0)*exp(-delta*z*FARADAY*v*(0.001)/R/(celsius+273)))
}
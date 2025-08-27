NEURON {
    SUFFIX sKCa
    USEION ca READ cai
    USEION k READ ki,ek WRITE ik
    RANGE  gk,isKCa
    GLOBAL sKCatau,rate_k,gmax_k
}

UNITS {
    (mM) = (milli/liter)
    (mA) = (milliamp)
    F = (faraday) (coulombs)	
}

PARAMETER {
    v (mV)
    dt (ms)
    gk = 0.0001 (mho/cm2)
    isKCa = 0.0 (mA/cm2)
    sKCatau = 6.1 (ms)
    ek 
    ki
    cai
    celsius	
}

ASSIGNED {
    ica (mA/cm2)
    ik (mA/cm2)
    winf 
    wtau (ms)
    rate_k
    gmax_k
}

STATE {
    w
}

BREAKPOINT {
    SOLVE integrate METHOD cnexp
    ik = (gk*gmax_k)*w*(v-ek)
    isKCa = ik
}

UNITSOFF

INITIAL {
    rate_k = 1.66
    gmax_k = 1.66
    setinf(cai)
    w = winf
}

DERIVATIVE integrate {
    setinf(cai)
    w' = (winf - w)/wtau
}

PROCEDURE setinf(cai) {
    LOCAL wcai
    
    wcai = cai*1000
    winf = 0.81/(1+exp((llog(wcai)+0.3)/-0.46))
    wtau = sKCatau/rate_k
}

FUNCTION llog(x) {  
    if (x>1e-11) {
        llog = log(x)
    }else{
        llog=0
    }
}

UNITSON
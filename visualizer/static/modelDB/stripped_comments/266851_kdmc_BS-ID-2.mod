NEURON {
    THREADSAFE
    SUFFIX kdmc
    USEION k READ ek WRITE ik
    RANGE  gbar, minf, mtau, hinf, htau, ik
    GLOBAL taumin
}

PARAMETER {
    gbar    = 0.002     (mho/cm2)

    celsius
    ek                  (mV)   
    v                   (mV)

    
    vhalfmt = -25       
    km      = 14        

    
    
    vhalfh  = -5        
    zetah   = 0.02      
    gmh     = 0.2       
    a0h     = 0.00058   
    taumin	= 0.1	(ms)		

    vhalfht = -100      
    kh      = 8         

    q10     = 3
}


UNITS {
    (mA) = (milliamp)
    (mV) = (millivolt)
    (pS) = (picosiemens)
    (um) = (micron)
}

ASSIGNED {
    ik      (mA/cm2)
    minf        mtau (ms)
    hinf        htau (ms)
}


STATE { m h }

BREAKPOINT {
    SOLVE states METHOD cnexp
    ik  = gbar*m*h*(v - ek)
}

INITIAL {
    trates(v)
    m   = minf
    h   = hinf
}

DERIVATIVE states {
    trates(v)
    m'  = (minf-m)/mtau
    h'  = (hinf-h)/htau
}

PROCEDURE trates(v) {
    LOCAL qt
    qt   = q10^((celsius-34)/10)

    minf = 1/(1 + exp(-(v-vhalfmt)/km))
    mtau = 1

    hinf = 1/(1 + exp((v-vhalfht)/kh))
    htau = exp(zetah*gmh*(v-vhalfh)) / (qt*a0h*(1 + exp(zetah*(v-vhalfh))))
    if(htau < taumin) { htau = taumin } 	
}
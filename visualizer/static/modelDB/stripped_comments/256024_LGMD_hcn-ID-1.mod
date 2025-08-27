UNITS {
    (mV) = (millivolt)
    (mA) = (milliamp)
    (S) = (siemens)
	(molar) = (1/liter)
	(mM) = (millimolar)
}



NEURON {
    THREADSAFE
    
    SUFFIX hcn

    NONSPECIFIC_CURRENT i
	USEION cn READ cni VALENCE 1

    
    
    
    RANGE gmax, tau, g, taumax, i

    
    GLOBAL e, taumin, vhalf, s1, s2
}

PARAMETER {
    
    gmax= 0.001 (S/cm2)
    e = -37 (mV)
    vhalf = -78.2 (mV)
    s1 = -12.5 (mV)
    s2 = 13 (mV)
    taumax = 1350 (ms)
    taumin = 10 (ms)
    
   	cnvm=15		(mV)
	lcp=2		(1)
	kD=6e-4		(mM)
	taucn=10	(ms)

}

ASSIGNED { 
    v (mV)
    cni (mM)
    
    i (mA/cm2)
   	vs	(mV)
    ninf
    tau (ms)
    g (S/cm2)
    
}

STATE {
    n
   	ov	(mV)
	ovs	(mV)
}

BREAKPOINT {
    SOLVE states METHOD cnexp
    n = 1/(1+(exp((vhalf-ov+ovs)/s1)))
    g  = gmax*n
    i  = gmax*n*(v-e)
}

INITIAL {
    settables(v)
    lci(cni)
    n = 1/(1+(exp((vhalf-v+vs)/s1)))
   	ov = v
	ovs = vs
}

DERIVATIVE states { 
    settables(v)
    lci(cni)
    
    ov' = (v-ov)/tau
	ovs' = (vs-ovs)/taucn
	
    
}


PROCEDURE settables(v (mV)) {
    
    

    TABLE tau DEPEND vhalf, s1, s2, taumax
          FROM -150 TO 50 WITH 750

    
    
    
    ninf = 1/(1+(exp((vhalf-v)/s1)))

    
    
    
    
    
    tau = 4*taumax/(1+exp((vhalf-v)/s2))*ninf+taumin

}


PROCEDURE lci(cni (mM)) { 
	TABLE vs DEPEND lcp, kD, cnvm
          FROM 0 TO 0.01 WITH 500
    
	vs = cnvm-cnvm/(1+(cni/kD)^lcp)
    
}
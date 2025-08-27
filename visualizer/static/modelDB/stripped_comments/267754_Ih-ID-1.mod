UNITS {
    (mA) = (milliamp)
    (mV) = (millivolt)
}
 
NEURON {
	THREADSAFE
    SUFFIX Ih
    USEION h READ eh WRITE ih VALENCE 1
    RANGE gkhbar,ih,g,t1,t2,t3,t4,t5,v_half,k
    GLOBAL rinf
}

INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}
 
PARAMETER {
    v (mV)
    p = 5 (degC)
    dt (ms)
	gkhbar = 0.000005 (mho/cm2)
    eh = -32.9 (mV)

	t1 = 1 (1/mV)
	t2 = -0.116 (1/mV)
	t3 = 1 (1/mV)
	t4 = 0.09 (1/mV)
	t5 = 100 (ms)
    
    v_half = -103.44 (mV)
    k = 8.63
}
 
STATE {
        r
}
 
ASSIGNED {
    ih (mA/cm2)
	rinf 
	tau_r (ms)
	g (siemens/cm2)
}
 
BREAKPOINT {
    SOLVE state METHOD cnexp
    ih = gkhbar*r*(v - eh)
	g = gkhbar*r
}
 
INITIAL {
	rates(v, v_half, k, t1, t2, t3, t4, t5)
	r = 0
}

DERIVATIVE state { 
	rates(v, v_half, k, t1, t2, t3, t4, t5)
	r' = (rinf - r)/tau_r
}
UNITSOFF

PROCEDURE rates(v(mV), v_half(mV), k, t1(1/mV), t2(1/mV), t3(1/mV), t4(1/mV), t5(ms)) {
    rinf = 1/(1 + exp((v-v_half)/k))
    tau_r = 1/(exp(-t1-t2*v) + exp(-t3+t4*v)) + t5
}
 
UNITSON
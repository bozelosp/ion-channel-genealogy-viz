UNITS {
        (mA) = (milliamp)
        (mV) = (millivolt)
}
 
NEURON {
        SUFFIX IMminret
        USEION k READ ek WRITE ik
        RANGE gbar,ik, vhalf1, vhalf2, k1, k2, c1, c2
        GLOBAL minf, mtau
}
 
INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}
 
PARAMETER {
        v (mV)
        dt (ms)
        gbar = 0.001 (mho/cm2)		
	ek
	mmin = 7	(ms)
	vhalf1 = -63
	vhalf2 = -63
	k1 = 15
	k2 = 15
	c1 = 0.003
	c2 = 0.003
}
 
STATE {
        m
}
 
ASSIGNED {
	ik
	minf 
	mtau	(ms)
}
 
BREAKPOINT {
        SOLVE state METHOD cnexp
        ik = gbar*m*(v - ek)
}
 
INITIAL {
	rates(v)
	m = minf
}

DERIVATIVE state { 
	rates(v)
	m' = (minf - m)/mtau
}

PROCEDURE rates(v(mV)) {  
                      

	LOCAL a, b
	minf = 1/(1 + exp(-(v+60(mV))/7(mV)))
	a = c1/exp(-(v-vhalf1)/k1)
	b = c2/exp((v-vhalf2)/k2)
	mtau = 1/(a+b)
	if (mtau<mmin) {mtau = mmin}
}
 
UNITSON
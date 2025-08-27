UNITS {
    (mV) = (millivolt)
    (mA) = (milliamp)
}

NEURON {
    THREADSAFE
    SUFFIX Kdr_dend
    USEION k READ ek WRITE ik
    RANGE gmax, gk
}

PARAMETER {
    gmax= 0.001 (mho/cm2)
    ek
    
	vhalf=-50
	tnmax=2.5
	tnmin=0.2
	tns=-30
	zn=15
}

ASSIGNED { 
    v (mV)
    ik (mA/cm2)
    ninf
    tau (ms)
    gk (S/cm2)
}

STATE {
    n
}

BREAKPOINT {
    SOLVE states METHOD cnexp
    gk  = gmax*n^4
    ik  = gk*(v-ek)
}

INITIAL {
    settables(v)
    n = ninf
}

DERIVATIVE states {  
    settables(v)      
   n' = (ninf - n)/tau
}

UNITSOFF

PROCEDURE settables(v (mV)) {

    TABLE ninf, tau DEPEND vhalf, zn, tns, tnmax
          FROM -150 TO 50 WITH 2000

    ninf = 1/(1+(exp((vhalf-v)/zn)))

    tau = tnmax/exp((vhalf-v)/tns)*ninf+tnmin
}

UNITSON
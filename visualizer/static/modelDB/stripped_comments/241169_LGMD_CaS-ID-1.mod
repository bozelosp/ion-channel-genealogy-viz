UNITS {
    (mV) = (millivolt)
    (mA) = (milliamp)
    (S) = (siemens)
}

NEURON {
    THREADSAFE
    SUFFIX CaS
    USEION ca READ eca WRITE ica
    RANGE gmax, g, taum, tauh
    GLOBAL vhalfm, tmax, vhalfh, hmin
}

PARAMETER {

    gmax = 0.001 (S/cm2)
    vhalfm =-30	(mV)
    vhalfh =-30	(mV)
    vn2= -40	(mV)
	s1 = 10.0	(mV)
	s2 = -8.0	(mV)
	tmax = 6	(ms)
	tmin = 2.5	(ms)
	kh = -5.5	(mV)
	hb = -130	(mV)
	hmin = 12	(ms)
	tadj = 1	(ms)
	hs = -30	(mV)
	nic = 0.0	(1)
	mp=2		(1)
}

ASSIGNED {
    v	(mV)
    eca (mV)
    
    minf (1)
    hinf (1)
    taum (ms)
    tauh (ms)
    ica	(mA/cm2)
	g	(S/cm2)
}

STATE {
    m
    h
}

INITIAL {
    settables(v)
    m = minf
    h = hinf
}

BREAKPOINT {
    SOLVE states METHOD cnexp
    g = gmax*m^mp*h
    ica  = g*(v-eca)
}


DERIVATIVE states { 
    settables(v)
	m' = (minf - m)/taum
	h' = (hinf - h)/tauh
}


PROCEDURE settables(v (mV)) {
    TABLE minf, taum, hinf, tauh DEPEND vhalfm, tmax, tmin, vhalfh, hmin, tadj
          FROM -100 TO 50 WITH 600

	
    minf = 1/(1 + exp((vhalfm-v)/s1))

    
	taum = tmin + 4*(tmax-tmin)/(1 + exp((vn2-v)/s2))*minf
	
	
	
	hinf = nic + (1-nic)/(1+exp((vhalfh-v)/kh))

    
	tauh = hmin + hinf*exp((hb-v)/hs)*tadj
	
}
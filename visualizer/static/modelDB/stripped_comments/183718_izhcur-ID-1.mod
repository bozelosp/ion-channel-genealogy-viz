NEURON {
	POINT_PROCESS izhcur
	NONSPECIFIC_CURRENT i_izh
	RANGE a,b,c,d,e,f,g,uinit,F,u,I
}

PARAMETER {
	a = 0.01	(1)
	b = 0.2		(1)
	c = -65		(mV)
	d = 2		(1)
	e = 0.04	(1)
	f = 5		(1)
	g = 140		(1)
	uinit = -14	(1)
	F = 1		(1) 
	I = 0		(mA/cm2) 
}

UNITS {
	(mV) = (millivolt)
	(mA) = (milliamp)
}

ASSIGNED {
	v (mV)
	i_izh (mA)
}

STATE { u }

BREAKPOINT {
	SOLVE states METHOD cnexp
	i_izh = (-1e-5)*F*(vinfi(v)-u)-I	
}

INITIAL {
	u = uinit
	net_send(0,1)					
}

DERIVATIVE states {
	UNITSOFF
	u'= F*a*(b*v-u)
	UNITSON
}

FUNCTION vinfi(v (mV)) {
	UNITSOFF
	vinfi = e*v*v + f*v + g
	UNITSON
}

NET_RECEIVE (w) {
	if (flag == 1) {
		WATCH (v > 30.0) 2
	} else {
		net_event(t)
		v = c
		u = u+d
  }
}
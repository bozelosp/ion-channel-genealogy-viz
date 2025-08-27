NEURON {
	POINT_PROCESS NMDA
	RANGE tau, time_interval, e, i,weight, NMDA_saturation_fact, flag, g
	NONSPECIFIC_CURRENT i
	GLOBAL gfac

	USEION nmda1 WRITE inmda1 VALENCE 0
	USEION nmda2 WRITE inmda2 VALENCE 0
	RANGE srcgid, targid, comp, synid
}

UNITS {
	(nA) = (nanoamp)
	(mV) = (millivolt)
	(uS) = (microsiemens)
	(mM) = (milli/liter)
}

PARAMETER {
	tau = 130.5 (ms)  <1e-9,1e9>	

	time_interval = 5 (ms) <1e-9,1e9>
	e = 0	(mV)
	weight = 2.5e-8 (uS)	
			 	
	NMDA_saturation_fact= 80e0 (1) 
		
		

	A_ = 0 (1) 
	BB1_ = 0 (1) 
	BB2_ = 0 (1) 
	Mg = 1.5 (mM) 
	gfac = 1
}

ASSIGNED {
	v (mV)
	i (nA)
	event_count (1)	
	k (uS/ms) 
	g (uS)
	A1_ (1)
	A2_ (1)
	B1_ (1)
	B2_ (1)
	Mg_unblocked (1)
	inmda1 (nA)
	inmda2 (nA)
	srcgid
	targid
	comp
	synid
}

STATE {
	A (uS)
	B (uS)
}

INITIAL {
	A_ =  exp(-2.847)  
	BB1_ = exp(-.693)  
	BB2_ = exp(-3.101) 
	g = 0
	A = 0
	B = 0
	k = 0
}

BREAKPOINT {
	SOLVE state METHOD cnexp
	g = A + B
	if (g > NMDA_saturation_fact * weight) { g = NMDA_saturation_fact * weight }
	g = g*Mg_unblocked*gfac
	i = g*(v - e)
	inmda1 = g
	inmda2 = -g
}

DERIVATIVE state {
	Mg_factor()
	B' = -B/tau
	A' = k
}

NET_RECEIVE(weight (uS)) {
	if (flag>=1) {
		
	
		k = k - weight/time_interval
	
		B = B + weight
		A = A - weight
	} else {
		
		net_send(time_interval, 1) 
	
		k = k + weight/time_interval


	}
}




PROCEDURE Mg_factor() {
UNITSOFF
           A1_ = exp(-.016*v - 2.91)
           A2_ = 1000.0 * Mg * exp (-.045 * v - 6.97)
           B1_ = exp(.009*v + 1.22)
           B2_ = exp(.017*v + 0.96)
UNITSON
           Mg_unblocked  = 1.0/(1.0 + (A1_+A2_)*(A1_*BB1_ + A2_*BB2_) /
                 (A_*A1_*(B1_+BB1_) + A_*A2_*(B2_+BB2_))  )
}
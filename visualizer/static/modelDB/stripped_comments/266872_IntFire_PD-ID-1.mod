NEURON {
	
	
	ARTIFICIAL_CELL IntFire_PD
	RANGE m, Iext, i
	
}

PARAMETER {
	
	taum = 10 (ms)
	
	tauref =2 (ms) 
	
	tausyn = 0.5 (ms) 
	
	Vreset = -65 (mV) 
	
	Vteta  = -50 (mV)
	
	Cm   = 250 (pF)
	Iext = 0(pA)
}


UNITS {
  (mV) = (millivolt)
  (pF) = (picofarad)
  (pA) = (picoamps)
}


ASSIGNED {
	m(mV)
	i(pA)
	t0(ms)
	refractory
}

FUNCTION M() {

}

FUNCTION I() {

}


INITIAL {
	t0 = t
	refractory = 0 
}

NET_RECEIVE (w) {
	if (refractory == 0) { 
		
		
		i=i-(i/tausyn)*(t-t0)
		i=i+w
		m=m+(((Vreset-m)/taum)+((i+Iext)/Cm))*(t-t0)
		
		
		
		
		
		
		
		
		
		
		
		t0 = t

		if (m > Vteta) {
			refractory = 1
			
			m = Vreset+30
			net_send(tauref, refractory)
			net_event(t)
		}
	}else if (flag == 1) { 
		t0 = t
		refractory = 0
		m = Vreset
	}
}
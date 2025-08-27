NEURON {
	POINT_PROCESS ScalExp2Syn
	
	
	RANGE tau1, tau2, Ee, i
	RANGE Enumsyns			
	RANGE id					
	RANGE Eintrvl, NEproc	
	RANGE Eg
	POINTER proc_num
	
	
	RANGE Egmax0, Egmax
	RANGE HSP_type			
	RANGE Vtrg, V
	RANGE Etau, Eenable
	RANGE continuous_update, event_window 
	RANGE order		
	POINTER Vsoma
	
	
	RANGE vavg, avgstrt
	RANGE mavg, mavgstrt, mavgintrvl

	NONSPECIFIC_CURRENT i
}

UNITS {
	(nA) = (nanoamp)
	(mV) = (millivolt)
	(uS) = (microsiemens)
	(nS) = (nanosiemens)
}

PARAMETER {
	tau1=.1 (ms) <1e-9,1e9>
	tau2 = 10 (ms) <1e-9,1e9>
	Ee=0	(mV)
	Enumsyns = 1
	id = 0
	Eintrvl = 500 (ms)		
	NEproc				
	
	
	Egmax0 = 1 (nS)
	Vtrg = -60 (mV)		
	Etau=1e5 (ms)			
	Eenable=1			
	order=0				

	continuous_update=1		
	event_window=5 (ms)		
	HSP_type = 1
	
	
	
	avgstrt=0 (ms)			
	mavgstrt=0 (ms)			
	mavgintrvl=5000 (ms)		
}

ASSIGNED {
	v (mV)
	i (nA)
	Eg (uS)
	factor
	Eorder

	vavg		  	
	stp1			
	stp2			
	mavg			
	vintrvl			
	n			

	in_window		
	synevent_time		
	proc_num
	Vsoma
	V
}

STATE {
	EA (uS)
	EB (uS)
	Egmax (nS)			
}

INITIAL {
	LOCAL j,tp

	if (tau1/tau2 > .9999) { tau1 = .9999*tau2 }
	EA = 0
	EB = 0
	tp = (tau1*tau2)/(tau2 - tau1) * log(tau2/tau1)
	factor = -exp(-tp/tau1) + exp(-tp/tau2)
	factor = 1/factor

	Egmax = Egmax0
	Eorder = Egmax^order
	if (HSP_type == 0) {
		V = Vsoma
	} else if (HSP_type == 1) {
		V = v
	}
	in_window = continuous_update	
	synevent_time = 0
	stp1 = 0
	stp2 = 0
	n = 1
	vavg = 0
	vintrvl = 0
	mavg = v
}

BREAKPOINT {
	Eg = EB - EA
	i = Eg*(v - Ee)
	SOLVE update METHOD after_cvode
	SOLVE state METHOD cnexp
}

DERIVATIVE state {
	EA' = -EA/tau1
	EB' = -EB/tau2

}

NET_RECEIVE(w) {
	LOCAL draw, j
	j = 0
	while (j < Enumsyns) {
		draw = floor(scop_random() * NEproc)	
		if (proc_num == draw) {
			synevent_time=t
			state_discontinuity(EA, EA + Egmax/1000*factor)
			state_discontinuity(EB, EB + Egmax/1000*factor)
			net_event(t)		
		}
		j = j + 1
	}
}

PROCEDURE update() {
	LOCAL sf
	
	
	if (Egmax<0) { Egmax=0 }
	Eorder = Egmax^order
	if (HSP_type == 0) {
		V = Vsoma
	} else if (HSP_type == 1) {
		V = v
	}
	sf = Vtrg - V
	if (Eenable) {

		Egmax = Egmax * (1 + exprand(sf) / Etau)
	}
	
	
	if (!continuous_update) {
		if (t-synevent_time<=event_window) { in_window=1
		} else { in_window=0 }
	}
	
	
	if (t>=avgstrt) {
   		vavg=(vavg*stp1+v)/(stp1+1)
		stp1=stp1+1
	}
	
	
	if (t>=mavgstrt) {
      		vintrvl=vintrvl+v
		stp2=stp2+1
   		if (t>mavgstrt+n*mavgintrvl-dt/2) { 
      			if (mavgintrvl>0) { mavg=vintrvl/stp2 }
      			n=n+1
      			vintrvl=0
			stp2=0
   		}
	}
}
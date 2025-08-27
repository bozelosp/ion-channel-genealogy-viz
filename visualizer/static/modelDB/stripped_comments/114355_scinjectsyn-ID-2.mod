NEURON {
	POINT_PROCESS ScalInjectSyn
	
	
	RANGE Ee, Eg, i
	
	
	RANGE Egmax0, Egmax
	RANGE HSP_type			
	RANGE Vtrg, V
	RANGE Etau, Eenable
	RANGE order		
	POINTER Vsoma
	
	
	RANGE vavg, avgstrt
	RANGE mavg, mavgstrt, mavgintrvl
	
	
	RANGE tau1, tau2, Enumsyns, id, Eintrvl, NEproc	
	POINTER proc_num
	RANGE continuous_update, event_window

	NONSPECIFIC_CURRENT i
}

UNITS {
	(nA) = (nanoamp)
	(mV) = (millivolt)
	(uS) = (microsiemens)
	(nS) = (nanosiemens)
}

PARAMETER {
	Ee=0	(mV)
	
	
	Egmax0 = 1 (nS)
	Vtrg = -60 (mV)		
	Etau=1e5 (ms)			
	Eenable=1			
	order=0				
	HSP_type = 1
	
	
	avgstrt=0 (ms)			
	mavgstrt=0 (ms)			
	mavgintrvl=5000 (ms)		
	
	
	tau1=.1 (ms) <1e-9,1e9>
	tau2 = 10 (ms) <1e-9,1e9>
	Enumsyns = 1
	id = 0
	Eintrvl = 500 (ms)
	NEproc
	continuous_update=1
	event_window=5 (ms)
}

ASSIGNED {
	v (mV)
	i (nA)
	Eg (uS)
	Eorder
	vavg		  	
	stp1			
	stp2			
	mavg			
	vintrvl			
	n			
	Vsoma
	V
	
	
	proc_num
}

STATE {
	Egmax (nS)			
}


INITIAL {
	Egmax = Egmax0
	Eorder = Egmax^order
	if (HSP_type == 0) {
		V = Vsoma
	} else if (HSP_type == 1) {
		V = v
	}
	stp1 = 0
	stp2 = 0
	n = 1
	vavg = 0
	vintrvl = 0
	mavg = v
}

BREAKPOINT {
	Eg = Egmax
	i = Eg*(v - Ee)
	SOLVE update METHOD after_cvode
	SOLVE state METHOD cnexp
}

DERIVATIVE state {
	Egmax' = Eenable * Eorder * (Vtrg - V) / Etau
}

NET_RECEIVE(w) {
	
}

PROCEDURE update() {
	
	
	if (Egmax<0) { Egmax=0 }
	Eorder = Egmax^order
	if (HSP_type == 0) {
		V = Vsoma
	} else if (HSP_type == 1) {
		V = v
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
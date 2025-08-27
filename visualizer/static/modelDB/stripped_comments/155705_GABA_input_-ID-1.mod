NEURON {
	POINT_PROCESS Gaba_synf
	RANGE e, i, g, area_cell 
	RANGE tau_d, frac_rec
	RANGE tau_1, tau_rec, tau_facil, U, u0
	NONSPECIFIC_CURRENT i
}

UNITS { 
	(nA)   = (nanoamp)
	(mV)   = (millivolt)
	(uS)   = (microsiemens)
	(um)   = (micron)
	}

PARAMETER {	
	
	
	tau_d = 10 (ms)      < 1e-9, 1e9 > 
	frac_rec = 0.9 		(1)  <0,1> 
	area_cell= 1 		(um2) 	   
	g=1 (1)
	
	

	
	
	e=-80	(mV)			 
	tau_1 = 1 (ms) < 1e-9, 1e9 >     
	tau_rec = 100 (ms) < 1e-9, 1e9 > 
	tau_facil =0
	u0 = 0 (1) < 0, 1 > 		 
	U = 0.04 (1) < 0, 1 > 		 
					 
	


}

ASSIGNED {
	v (mV)
	i (nA)
	x
	
}

STATE {
geff (siemens)	
}

INITIAL {	
geff=0	 
}

BREAKPOINT {
	SOLVE state METHOD cnexp
	i=geff*(v - e)

}

DERIVATIVE state {
	geff'=-geff/tau_1
}

NET_RECEIVE(w (us),y, z, u,tp(ms)) {
INITIAL {

	y = 0
	z = 0
	u = u0 
	tp = t
}
	
	
	z = z*exp(-(t - tp)/tau_rec)
	z = z + ( y*(exp(-(t - tp)/tau_1) - exp(-(t - tp)/tau_rec)) / ((tau_1/tau_rec)-1) )
	
	y = y*exp(-(t - tp)/tau_1)
	x = 1-y-z
	

	
	if (tau_facil > 0) {u = u*exp(-(t - tp)/tau_facil)} else {u = U}
	if (tau_facil > 0) {u= u + U*(1-u)}			
	geff=geff + w*x
	y= y + x
	tp = t






}
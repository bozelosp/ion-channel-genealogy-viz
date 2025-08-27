NEURON {
	POINT_PROCESS Gaba_syn
	RANGE e, i, g
	RANGE tau_d, frac_rec
	RANGE area_cell
	NONSPECIFIC_CURRENT i
}

UNITS {
	(nA)   = (nanoamp)
	(mV)   = (millivolt)
	(uS)   = (microsiemens)
	(um)   = (micron)
}

PARAMETER {
	
	
	
	e=-80	(mV)
	
	
		tau_d = 10 (ms) 	< 1e-9, 1e9 > 
	    frac_rec = 0.9 (1) 	<0,1>
	
	
	   g=1e-5 (uS/um2) 
                   
                   
                   
    
        area_cell= 1 (um2)
}

ASSIGNED {
	v (mV)
	i (nA)
	g_eff (uS)
	
}

STATE {
	s
}

INITIAL {
	s=0 
	g_eff=g*area_cell
}

BREAKPOINT {
	SOLVE state METHOD cnexp
	i =g_eff*s*(v - e) 
}

DERIVATIVE state {
	s' = 	-s/tau_d
}

NET_RECEIVE(weight,s_tp,tp(ms)) {
	
	UNITSOFF
	
    
    s_tp =s_tp*weight*exp(-(t-tp)/tau_d) 
    
	UNITSON
	
	
	
	 s=s + frac_rec*weight*(1-s_tp)	
    s_tp = s_tp + frac_rec*(1-s_tp)
	
    tp=t   
}
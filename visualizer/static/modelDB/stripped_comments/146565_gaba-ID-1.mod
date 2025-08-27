NEURON {
	POINT_PROCESS gaba
	
	NONSPECIFIC_CURRENT i

	RANGE e ,gmax,local_v,i
	RANGE del,Tspike,Nspike
	RANGE taudgaba,dgaba,decaygaba
	RANGE R,D
	GLOBAL risetime,decaytime
}

PARAMETER {
	gmax=.5	(nS)
	e= -70.0	(mV)
	risetime=1	(ms)	
	decaytime=20(ms)	

	v		(mV)
	del=30	(ms)
	Tspike=10	(ms)
	Nspike=1
	taudgaba=200	(ms)
	decaygaba=0.6
	}

ASSIGNED { 
	i		(nA)  
	local_v	(mV)
}
STATE {
	dgaba
	R
	D
}

INITIAL {
      dgaba=1 
	R=0
	D=0
}    

BREAKPOINT {  
    
	LOCAL count
	SOLVE state METHOD cnexp
	FROM count=0 TO Nspike-1 {
		IF(at_time(count*Tspike+del)){
			state_discontinuity( R, R+ gmax*dgaba)
			state_discontinuity( D, D+ gmax*dgaba)
			state_discontinuity( dgaba, dgaba* decaygaba)
		}
	}
	i= (D-R)*(1e-3)* (v- e)

	
	local_v=v
}

DERIVATIVE state {
	R'=-R/risetime
	D'=-D/decaytime
	dgaba'=(1-dgaba)/taudgaba

}
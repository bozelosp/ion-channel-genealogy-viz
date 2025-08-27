NEURON {
	POINT_PROCESS gaba_syn
	NONSPECIFIC_CURRENT i
	RANGE i,ggaba
	RANGE decaygaba,dgaba,taudgaba
	RANGE facilgaba,fgaba,taufgaba
	
	RANGE risetime,decaytime,e 
}
PARAMETER {
	e= -60.0	(mV)
	risetime=1	(ms)	
	decaytime=20(ms)	

	v		(mV)
	taudgaba=200	(ms)
	decaygaba=0.8
	taufgaba=200    (ms)
	facilgaba=0.0
}
ASSIGNED {
	i		(nA)  
	ggaba
    factor     
}

STATE {
	dgaba
	fgaba
	R
	D
}

INITIAL {
	LOCAL tp
    dgaba=1 
    fgaba=1
	R=0
	D=0
	ggaba=0
	
    tp = (risetime*decaytime)/(decaytime - risetime) * log(decaytime/risetime)
    factor = -exp(-tp/risetime) + exp(-tp/decaytime)
    factor = 1/factor
}
BREAKPOINT {
	SOLVE state METHOD cnexp
	ggaba=D-R
	i=(1e-3)*ggaba*(v-e)
}
NET_RECEIVE(weight) {
    R = R + factor*weight*dgaba*fgaba
    D = D + factor*weight*dgaba*fgaba
    dgaba = dgaba* decaygaba
    fgaba = fgaba + facilgaba
}
DERIVATIVE state {
	R'=-R/risetime
	D'=-D/decaytime
	dgaba'=(1-dgaba)/taudgaba
	fgaba'=(1-fgaba)/taufgaba

}
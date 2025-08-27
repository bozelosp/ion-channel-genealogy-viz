NEURON {
	SUFFIX iapnew
	USEION na READ ena WRITE ina
	USEION k READ ek WRITE ik
	GLOBAL inf
	RANGE gnabar, gkbar, ena, ek, ik
}

UNITS {
	(mA) = (milliamp)
	(mV) = (millivolt)
}

INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}
PARAMETER {
	v (mV)
	celsius = 37	(degC)
	gnabar=0.0 (mho/cm2)
	gkbar=.12 (mho/cm2)
	
	
	naactvha = 40 (mV)
	naiactvha = 45 (mV)
	kactvha = 40 (mV)	
	naactslope = -3 (mV)
	nainactslope = 3 (mV)
	kactslope = -3 (mV)
	Nainactivationtau = 0.5
}

STATE {
	m h n
}

ASSIGNED {
	ena (mV)
	ek (mV)
	ina (mA/cm2)
	ik (mA/cm2)
	inf[3]
	tau[3]
}


BREAKPOINT {
	SOLVE states METHOD cnexp
	ina = gnabar*m*m*h*(v - ena)
	ik = gkbar*n*n*(v - ek)
}

DERIVATIVE states {	
	mhn(v*1(/mV))
	m' = (inf[0] - m)/tau[0]
	h' = (inf[1] - h)/tau[1]
	n' = (inf[2] - n)/tau[2]
}

FUNCTION varss(v, i) {
	if (i==0) {
		varss = 1 / (1 + exp((v + naactvha)/(naactslope))) 
	}
	else if (i==1) {
		varss = 1 / (1 + exp((v + naiactvha)/(nainactslope))) 
	}
	else {
		
		varss = 1 / (1 + exp((v + kactvha)/(kactslope))) 
	}
}

FUNCTION vartau(i) {
	if (i==0) {
		vartau = 0.05  
	}
	else if (i==1) {
		vartau = Nainactivationtau   
	}
	else {
		vartau = 2     
	}
}

PROCEDURE mhn(v) {LOCAL a, b 
	TABLE inf,tau 
	DEPEND celsius, dt, kactslope, kactvha, nainactslope,naiactvha,naactslope, naactvha, Nainactivationtau
	FROM -100 TO 100 WITH 2000 
	FROM i=0 TO 2 {
		tau[i] = vartau(i)
		inf[i] = varss(v,i)
	}
}
UNITSON
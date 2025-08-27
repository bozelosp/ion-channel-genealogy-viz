NEURON 
{
	SUFFIX stp
   	USEION ca READ cai
	RANGE TTINF, tauTT, fT, tauTC
}

UNITS 
{
	(nA) = (nanoamp)
	(mV) = (millivolt)
	(uS) = (microsiemens)
	(molar) = (1/liter)
	(mM) = (millimolar)
	R = (k-mole) (joule/degC)
	FARADAY = (faraday) (coulomb)
    	(mA) = (milliamp)
    	(S) = (siemens)
}

PARAMETER 
{
	TTINF=20 	(mM) <1e-9,1e9> 
	tauTT=2000   	(ms) <1e-9,1e9> 
	
	fT=2000		(/mM/mM/mM/ms)	

	
	tauTC=0.001	(ms)
}


ASSIGNED 
{
   	cai     (mM)
}

STATE
{
	TT 	(mM) 
	TC	(mM) 
}

INITIAL 
{
	TT = TTINF 
	TC = 0	   
}

BREAKPOINT 
{
	SOLVE states METHOD cnexp

}

DERIVATIVE states
{
	TT'=(TTINF-TT)/tauTT - fT*TT*(cai-50e-6)*(cai-50e-6)*(cai-50e-6)
	TC'=fT*TT*(cai-50e-6)*(cai-50e-6)*(cai-50e-6) - TC/tauTC
}
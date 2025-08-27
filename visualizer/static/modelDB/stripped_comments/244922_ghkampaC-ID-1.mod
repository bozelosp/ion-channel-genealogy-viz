NEURON 
{
	POINT_PROCESS ghkampaC
	USEION na WRITE ina
	USEION k WRITE ik
		
	RANGE C, TRise, tau, lr
	RANGE Alpha, Beta

	RANGE iampa
	RANGE P, Pmax
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
	
}

PARAMETER 
{
	TRise=0.6 	(ms)<1e-9,1e9>  
	tau=3   	(ms)<1e-9,1e9>	 
	
	nai = 18	(mM)	
	nao = 140	(mM)
	ki = 140	(mM)	
	ko = 5		(mM)
	celsius	 = 35	(degC)
	Pmax = 1e-6
}


ASSIGNED 
{
	v (mV)
	ina     (nA)
	ik      (nA)
	Alpha	(/ms mM)
	Beta	(/ms)	
	C		(mM)		
	P (cm/s)
	iampa	(nA)
	Area (cm2)
	lr

}

STATE
{
	S				
}

INITIAL 
{
	S = 0
	Beta=1/tau
	Alpha=1/TRise - Beta
	Area=1
	
}

BREAKPOINT 
{
	SOLVE states METHOD cnexp
	P = (Pmax * S * (Alpha+Beta)) / (Alpha*(1-1/exp(1)))

	ina= P*ghk(v, nai, nao,1) * Area
	ik= P*ghk(v, ki, ko,1)* Area
	iampa=ina+ik 
}

DERIVATIVE states
{	
	S'=Alpha * C * (1-S) - Beta * S
}

FUNCTION ghk(v(mV), ci(mM), co(mM),z) (0.001 coul/cm3) {
	LOCAL arg, eci, eco
	arg = (0.001)*z*FARADAY*v/(R*(celsius+273.15))
	eco = co*efun(arg)
	eci = ci*efun(-arg)
	ghk = (0.001)*z*FARADAY*(eci - eco)
}

FUNCTION efun(z) {
	if (fabs(z) < 1e-4) {
		efun = 1 - z/2
	}else{
		efun = z/(exp(z) - 1)
	}
}
NEURON {
	POINT_PROCESS Wghkampa
	USEION na WRITE ina
	USEION k WRITE ik
	USEION ca READ cai	
	
	RANGE taur, taud
	RANGE iampa,winit
	RANGE P, Pmax, lr
}

UNITS {
	(nA) = (nanoamp)
	(mV) = (millivolt)
	(uS) = (microsiemens)
	(molar) = (1/liter)
	(mM) = (millimolar)
	FARADAY = (faraday) (coulomb)
	R = (k-mole) (joule/degC)
  	(um)    = (micron)
  	PI      = (pi)       (1)

}

PARAMETER {
	taur=2 		(ms) <1e-9,1e9>
	taud = 10 	(ms) <1e-9,1e9>
	nai = 18	(mM)	
	nao = 140	(mM)
	ki = 140	(mM)	
	ko = 5		(mM)
	celsius		(degC)
	Pmax=1e-6   (cm/s)	
	alpha1=0.35	
	beta1=80
	alpha2=0.55
	beta2=80
	winit=1		(1)
}

ASSIGNED {
	ina     (mA/cm2)
	ik      (mA/cm2)
	v (mV)
	P (cm/s)
	factor
	iampa	(mA/cm2)
	lr
	cai			(mM)
	Area (cm2)
	diam (um)
}

STATE {
	A (cm/s)
	B (cm/s)
	w (1)
}

INITIAL {
	LOCAL tp
	if (taur/taud > .9999) {
		taur = .9999*taud
	}
	A = 0
	B = 0
	tp = (taur*taud)/(taud - taur) * log(taud/taur)
	factor = -exp(-tp/taur) + exp(-tp/taud)
	factor = 1/factor
	Area=PI*diam*1e-2 
	w=winit
}

BREAKPOINT {
	SOLVE state METHOD cnexp
	P=B-A



	ina = P*w*ghk(v, nai, nao,1)/Area	
	ik = P*w*ghk(v, ki, ko,1)/Area
	iampa = ik + ina
}

DERIVATIVE state {
	lr=eta(cai)
	w' = lr*(Omega(cai)-w)
	A' = -A/taur
	B' = -B/taud
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

FUNCTION eta(ci (mM)) { 
	LOCAL inv, P1, P2, P3, P4
	P1=100	
	P2=P1*1e-4	
	P4=1e3
	P3=3		

	ci=(ci-1e-4)*1e3 	

	inv=P4 + P1/(P2+ci*ci*ci) 
	eta=1/inv
}	

FUNCTION Omega(ci (mM)) {
	ci=(ci-1e-4)*1e3	
	Omega=0.25+1/(1+exp(-(ci-alpha2)*beta2))-0.25/(1+exp(-(ci-alpha1)*beta1))
}
	
NET_RECEIVE(weight (uS)) { 	
							
							
	state_discontinuity(A, A + Pmax*factor)
	state_discontinuity(B, B + Pmax*factor)
}
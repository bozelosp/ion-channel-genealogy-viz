NEURON {
	POINT_PROCESS ghknmda
	USEION na WRITE ina
	USEION k WRITE ik
	USEION ca READ cai, cao WRITE ica
	
	RANGE taur, taud
	RANGE inmda

	RANGE P, mg, Pmax 
	GLOBAL  mgb
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
	taur=5 (ms) <1e-9,1e9>
	taud = 50 (ms) <1e-9,1e9>
	cai = 100e-6(mM)	
	cao = 2		(mM)
	nai = 18	(mM)	
	nao = 140	(mM)
	ki = 140	(mM)	
	ko = 5		(mM)
	celsius		(degC)
	mg = 2		(mM)
	Pmax=1e-6   (cm/s)	
						
	
}

ASSIGNED {
	ina     (mA/cm2)
	ik      (mA/cm2)
	ica     (mA/cm2)
	v (mV)
	P (cm/s)
	factor
	mgb	(1)
	inmda	(mA/cm2)

	diam (um)
}

STATE {
	A (cm/s)
	B (cm/s)
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
}

BREAKPOINT {
	SOLVE state METHOD cnexp
	P=B-A
	mgb = mgblock(v)


 
	ina = P*mgb*ghk(v, nai, nao,1)     
	ica = P*10.6*mgb*ghk(v, cai, cao,2)  
	ik = P*mgb*ghk(v, ki, ko,1)      
	inmda = ica + ik + ina
}

DERIVATIVE state {
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

FUNCTION mgblock(v(mV)) (1){
	TABLE 
	DEPEND mg
	FROM -140 TO 80 WITH 1000 

	

	mgblock = 1 / (1 + exp(0.062 (/mV) * -v) * (mg / 3.57 (mM)))
}

NET_RECEIVE(weight (uS)) { 	
							
							
	state_discontinuity(A, A + Pmax*factor)
	state_discontinuity(B, B + Pmax*factor)
}
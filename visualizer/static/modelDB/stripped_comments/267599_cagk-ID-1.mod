UNITS {
	(molar) = (1/liter)
}

UNITS {
	(mV) =	(millivolt)
	(mA) =	(milliamp)
	(mM) =	(millimolar)
}


NEURON {
	SUFFIX mykca
	USEION k READ ek WRITE ik
	USEION ca READ cai
	RANGE gkbar,ik, cainit
	GLOBAL oinf, tau
	
}

UNITS {
	FARADAY = (faraday)  (kilocoulombs)
	R = 8.313424 (joule/degC)
}

PARAMETER {
      v		(mV)
	dt		(ms)
	ek		(mV)
	cai (mM)
	celsius = 20	(degC)
	gkbar = 0.01	(mho/cm2)	
	cainit = 100e-6 (mM)

      d1 =1
     	d2 = 1.5
	k1 = 0.18	(mM)
	k2 = 0.011	(mM)
	bbar = 0.28	(/ms)
	abar = 0.48	(/ms)


	
	
	
	
	
	


        st=1            (1)
}

ASSIGNED {
	ik		(mA/cm2)
	oinf
	tau		(ms)
	
      
}

INITIAL {
        rate(v,cainit)
        o=oinf
}

STATE {	o }		

BREAKPOINT {
	SOLVE state METHOD cnexp
	ik = gkbar*o^st*(v - ek)
}

DERIVATIVE state {	
	rate(v, cai)
	o' = (oinf - o)/tau
}

FUNCTION MyExp(x) {
    if (x<-50) {MyExp=0}
    else if (x>50) {MyExp=exp(50)}
    else {MyExp=exp(x)}
}

FUNCTION alp(v (mV), c (mM)) (1/ms) { 
	alp = c*abar/(c + exp1(k1,d1,v))
}

FUNCTION bet(v (mV), c (mM)) (1/ms) { 
	bet = bbar/(1 + c/exp1(k2,d2,v))
}

FUNCTION exp1(k (mM), d, v (mV)) (mM) { 
	exp1 = k*MyExp(-2*d*FARADAY*v/R/(273.15 + celsius))
}

PROCEDURE rate(v (mV), c (mM)) { 
	LOCAL a
	a = alp(v,c)
	tau = 1/(a + bet(v, c))
	oinf = a*tau
	
}
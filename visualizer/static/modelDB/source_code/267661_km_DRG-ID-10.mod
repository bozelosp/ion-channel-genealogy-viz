TITLE DRG KM channel based on Barkai et al. 2017 biophysical properties
: adaptation of M. Migliore June 2006 

UNITS {
	(mA) = (milliamp)
	(mV) = (millivolt)

}

NEURON {
	SUFFIX km
	USEION k READ ek WRITE ik
    RANGE  gbar, ik, vhalfl, localtemp
    GLOBAL inf, tau
}

PARAMETER {
	v 		(mV)
	ek
	localtemp = 37
	gbar=.0001 	(mho/cm2)
        vhalfl=-42 :40   	(mV)
		kl=12 :-6 :-10
        vhalft=-42   	(mV)
        a0t=0.009      	(/ms)
        zetat=7    	(1)
        gmt=.4   	(1)
	b0=60
	st=1
	A1= -0.07245
    A2= 1.13462
}

STATE {
        m
}

ASSIGNED {
	ik (mA/cm2)
        inf
	tau
        taua
	taub
	q10
}

INITIAL {
	rate(v)
	m=inf
	q10=5^((localtemp - 37)/10)
}


BREAKPOINT {
	SOLVE state METHOD cnexp
	ik = gbar*m^st*(v-ek)
}


FUNCTION alpt(v(mV)) {
  alpt = exp(0.0378*zetat*(v-vhalft)) 
}

FUNCTION bett(v(mV)) {
  bett = exp(0.0378*zetat*gmt*(v-vhalft)) 
}

DERIVATIVE state {
        rate(v)
:        if (m<inf) {tau=taua} else {tau=taub}
	m' = q10 * (inf - m)/tau
}

PROCEDURE rate(v (mV)) { :callable from hoc
        LOCAL a
        :inf = (1/(1 + exp((v-vhalfl)/kl)))
        
		inf = A2 + (A1-A2)*(1/(1 + exp((v-vhalfl)/kl)))
		
		a = alpt(v)
  :tau=110+v
        tau = b0 + bett(v)/(a0t*(1+a))
        taua = 50
        taub = 300
}

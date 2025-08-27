TITLE CA1 KM channel from Mala Shah
: M. Migliore June 2006

UNITS {
	(mA) = (milliamp)
	(mV) = (millivolt)

}

PARAMETER {
	v 		(mV)
	ek
	celsius 	(degC)
	gbar=.0001 	(mho/cm2)
        vhalfl=-40   	(mV)
	kl=-10
        vhalft=-42   	(mV)
        a0t=0.003      	(/ms)
        zetat=7    	(1)
        gmt=.4   	(1)
        vhalftb=-42   	(mV)
        a0tb=0.003      	(/ms)
        zetatb=7    	(1)
        gmtb=.4   	(1)
	q10=5
	b0=60
	b0b=2000
	st=1
	sh =0
}


NEURON {
	SUFFIX kmb
	USEION k READ ek WRITE ik
        RANGE  gbar,ik, sh
      GLOBAL inf, tau, taua, taub
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
}

INITIAL {
	rate(v)
	m=inf
}


BREAKPOINT {
	SOLVE state METHOD cnexp
	ik = gbar*m^st*(v-ek)
}


FUNCTION alpt(v(mV)) {
  alpt = exp(0.0378*zetat*(v-vhalft-sh)) 
}

FUNCTION bett(v(mV)) {
  bett = exp(0.0378*zetat*gmt*(v-vhalft-sh)) 
}

FUNCTION alptb(v(mV)) {
  alptb = exp(0.0378*zetatb*(v-vhalftb-sh)) 
}

FUNCTION bettb(v(mV)) {
  bettb = exp(0.0378*zetatb*gmtb*(v-vhalftb-sh)) 
}

DERIVATIVE state {
        rate(v)
    if (m<inf) {tau=taua} else {tau=taub}
	m' = (inf - m)/tau
}

PROCEDURE rate(v (mV)) { :callable from hoc
        LOCAL a,qt
        qt=q10^((celsius-35)/10)
        inf = (1/(1 + exp((v-vhalfl-sh)/kl)))
        taua = (b0 + bett(v)/(a0t*(1+alpt(v))))/qt
        taub = (b0b + bettb(v)/(a0tb*(1+alptb(v))))/qt

}















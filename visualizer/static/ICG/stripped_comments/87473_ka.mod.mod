UNITS {
	(mA) = (milliamp)
	(mV) = (millivolt)

}

PARAMETER {
	celsius 	(degC)
	gbar=.004 (mho/cm2)
	vha=-40 (mV)
	vhb=-70 (mV)
	aa=0.05	(/mV)
	ab=-0.1	(/mV)
	btau=10 (ms)
        v       (mV)
        ek      (mV)
	basic = 0
}


NEURON {
	SUFFIX ka
	USEION k READ ek WRITE ik
        RANGE gbar,gka
        RANGE ainf, binf, btau
	RANGE tot
}

STATE {
	b
}

ASSIGNED {
	ik (mA/cm2)
	tot (mA/cm2)
        gka  (mho/cm2)
	ainf
        binf
}

INITIAL {
        rates(v)
        b=binf

}

BREAKPOINT {
	SOLVE state METHOD cnexp
	gka = gbar*ainf*b
	tot = gka*(v-ek)
	ik = gka*(v-ek)

}

FUNCTION expn(v (mV),a(/mV), vhalf(mV)) {
  	expn = exp(-2*a*(v-vhalf))
}

DERIVATIVE state {     
        rates(v)
        b' = (binf - b)/btau
}

PROCEDURE rates(v (mV)) { 
	binf = 1/(1 + expn(v,ab,vhb))
	ainf = 1/(1 + expn(v,aa,vha))
	if( basic > 0 ) {
		
	 	ainf = 1
	}
}
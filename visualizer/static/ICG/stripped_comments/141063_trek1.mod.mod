NEURON {
	SUFFIX trek1
	USEION k READ ek WRITE ik
	RANGE gkbar
}


UNITS {
	(mA) = (milliamp)
	(mV) = (millivolt)
}

PARAMETER {
	gkbar   = .005  (mho/cm2)
	Pmin = 0.1  (1)    
	q10 = 5   (1)
	Tm  =  39   (degC)
}   

ASSIGNED {
	ek      (mV)
	celsius	(degC)
	ik      (mA/cm2)
	v	(mV)
	o
}

BREAKPOINT {
	Ps()
    ik  = o * gkbar * (v - ek)
}

INITIAL {
    o = Pmin + (1-Pmin)/(1+exp(-(celsius-Tm)*q10/10 (degC))) 
	ik  = o * gkbar * (v - ek)
}

PROCEDURE Ps(){
    o = Pmin + (1-Pmin)/(1+exp(-(celsius - Tm)*q10/10 (degC))) 
}
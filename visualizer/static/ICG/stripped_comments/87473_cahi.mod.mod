UNITS {
	(molar) =	(1/liter)
	(mM) =	(millimolar)
	(mA) = (milliamp)
	(mV) = (millivolt)

}

NEURON {
	SUFFIX cahi
	USEION ca READ eca WRITE ica
        RANGE gbar
        GLOBAL xinf
	RANGE tot
}

PARAMETER {
	v (mV)
	celsius 	(degC)
	
	gbar = 0.0002 (mho/cm2)
	
	
	xtau=5 (ms)
	Kc=1   (mM)
	ax=0.08 (/mV)
	vhx=-30 (mV)
	vrest = 124	(mV)
	simp = 0
}


STATE {
	x
}

ASSIGNED {
	ica (mA/cm2)
	tot (mA/cm2)
	cai (mM)
        gca (mho/cm2)
        xinf
	eca (mV)
}

INITIAL {
	rate(v)
	x = xinf

}

BREAKPOINT {
	SOLVE state METHOD cnexp
	gca = gbar*x*x*Kc/(Kc+cai)
	tot = gca*(v-eca)
	if( simp > 0 ) {
		
		gca = gbar*x
		
	}
	
	ica = gca*(v-eca)
}

FUNCTION expn(v (mV),a(/mV), vhalf(mV)) {
  	expn = exp(-2*a*(v-vhalf))
}

DERIVATIVE state {  
        rate(v)
        x' = (xinf - x)/xtau
}

PROCEDURE rate(v (mV)) { 
	xinf = 1/(1 + expn(v,ax,vhx))
}
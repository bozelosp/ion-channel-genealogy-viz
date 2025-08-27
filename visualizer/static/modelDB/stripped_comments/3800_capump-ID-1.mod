NEURON {
	SUFFIX ca_pump
	USEION ca READ cao, cai WRITE ica	
	RANGE IpCamax, ica, ina, ipca
	GLOBAL dummy 
}

UNITS {
	(mA) = (milliamp)
	(mV) = (millivolt)
}

PARAMETER {
	IpCamax = 0.55e-3       (mA/cm2) <0,1e6>
	
}

ASSIGNED {
	celsius (degC)
	v (mV)
	ica (mA/cm2)
	cao (mM)
        cai (mM)
	ipca (mA/cm2) 
	dummy
}


BREAKPOINT { LOCAL q10
	
		q10 = 3^((celsius - 37)/10 (degC))
	
		ipca = q10*IpCamax*cai/(0.0005 + cai)
	
		ica =  ipca
}
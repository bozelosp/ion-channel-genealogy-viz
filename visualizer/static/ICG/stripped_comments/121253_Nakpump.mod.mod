NEURON {
	SUFFIX NaKpump
	USEION k READ ko WRITE ik
	USEION na WRITE ina
	RANGE ik, ina , INaKmax, ink, Kmko, Kmnai
	GLOBAL dummy 
}

UNITS {
	(mA) = (milliamp)
	(mV) = (millivolt)
	
}

PARAMETER {
	INaKmax = 2.45e-3 (mA/cm2) <0,1e6>
	Kmnai = 27.9 (mM)    <0,1e6>
	Kmko = 5.3 (mM)    <0,1e6>
	nai = 10 (mM)
}

ASSIGNED {
	celsius (degC)
	v (mV)
	ko (mM)
	ik (mA/cm2)
	ina (mA/cm2)
	ink (mA/cm2)
	dummy
}

BREAKPOINT { LOCAL q10 , fnk 
	
	q10 = 3^((celsius - 37)/10 (degC))
	
	fnk = (v + 150)/(v + 200)
				
	ink= q10*INaKmax*fnk/((1 + (Kmnai/nai)^1.5)*(1 + Kmko/ko))
	ina = 3*ink
	ik = -2*ink
}
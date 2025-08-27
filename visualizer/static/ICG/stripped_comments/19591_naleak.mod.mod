UNITS {
	(mV) = (millivolt)
	(mA) = (milliamp)
}

INDEPENDENT { v FROM -100 TO 50 WITH 50	(mV) }

NEURON {
	SUFFIX naleak
	
	USEION na READ ena WRITE ina
        RANGE gbar
}

PARAMETER {
	gbar = 1.167e-5	(mho/cm2)
	
}

ASSIGNED {
    ina    (mA/cm2)
    ena (mV)
}

BREAKPOINT {
	ina = gbar*(v - ena)
}
UNITS {
	(mV) = (millivolt)
	(mA) = (milliamp)
}

INDEPENDENT { v FROM -100 TO 50 WITH 50	(mV) }

NEURON {
	SUFFIX caleak
	USEION ca READ eca WRITE ica
	RANGE gbar, i
        
}

PARAMETER {
	gbar = 1.616e-6	(mho/cm2)
	
}

ASSIGNED { 
           eca (mV)
           i    (mA/cm2)
           ica	(mA/cm2)}

BREAKPOINT {
	i = gbar*(v - eca)
        ica = i
}
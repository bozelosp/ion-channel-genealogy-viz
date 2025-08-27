NEURON {
        SUFFIX NaP
        USEION na READ ena WRITE ina 
        RANGE g, ina
}

UNITS {
	(S)  = (siemens)
        (mA) = (milliamp)
        (mV) = (millivolt)
}

PARAMETER {
        g	(S/cm2)
}

ASSIGNED {
        v	(mV)
	ena	(mV)
        ina	(mA/cm2)
	minf
}

BREAKPOINT { 
	rates()
	ina= g* minf* (v- ena) 
}

INITIAL {
	rates()
}

PROCEDURE rates() { UNITSOFF
	minf= 1/ (1+ exp(-(v+ 40)/ 5))
} UNITSON
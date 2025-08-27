UNITS {
        (mA) = (milliamp)
        (mV) = (millivolt)
}

NEURON {
        SUFFIX ICa
        USEION ca READ eca WRITE ica
        RANGE gcabar, gca, m
}
 
PARAMETER {
        gcabar 	= 0.0015 (mho/cm2)	<0,1e9>	
	
}
 
ASSIGNED {
        eca (mV)
        v 	(mV)
	gca 	(mho/cm2)
        ica 	(mA/cm2)
	m	(1)
}
 
BREAKPOINT {
	UNITSOFF
        m =  1 / ( 1 +  exp(-(v+20)/9) )
	UNITSON
        gca = gcabar*m*m
	ica = gca * ( v - eca )
}
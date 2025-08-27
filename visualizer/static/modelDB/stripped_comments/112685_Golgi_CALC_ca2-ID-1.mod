NEURON {
        SUFFIX Golgi_CALC_ca2
        USEION ca2 READ ica2, ca2o WRITE ca2i VALENCE 2
        RANGE d, beta, ca2i0, ca2_pump_i
}

UNITS {
        (mV)    = (millivolt)
        (mA)    = (milliamp)
	(um)    = (micron)
	(molar) = (1/liter)
        (mM)    = (millimolar)
   	F      = (faraday) (coulomb)
}

PARAMETER {
        ica2             (mA/cm2)
        celsius    (degC)
        d = .2          (um)
        ca2o        (mM)         
        ca2i0     (mM)         
        beta = 1.3        (/ms)
}
ASSIGNED {
	ca2_pump_i	(mA)
}

STATE {
	ca2i (mM)
}

INITIAL {
        ca2i = ca2i0 
}

BREAKPOINT {
       SOLVE conc METHOD derivimplicit
}

DERIVATIVE conc {    
	
	ca2_pump_i = 2*beta*(ca2i-ca2i0)
	
	ca2i' = -ica2/(2*F*d)*(1e4) - beta*(ca2i-ca2i0)
}
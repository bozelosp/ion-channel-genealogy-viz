NEURON {
        SUFFIX Golgi_CALC_ca2
        USEION ca2 READ ica2, ca2o WRITE ca2i VALENCE 2
	RANGE Q10_diff,beta_Q10, ic, fix_celsius
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
        fix_celsius = 37 (degC)
        d = .2          (um)
        ca2o        (mM)
        ca2i0     (mM)
        beta = 1.3        (/ms)
	Q10_diff = 1.7
}
ASSIGNED {
	beta_Q10 (mho/cm2)
	ca2_pump_i	(mA)
	ic
  tau (ms)
}

STATE {
	ca2i (mM)
}

INITIAL {
	beta_Q10 = beta*(Q10_diff^((fix_celsius-23)/10))
        ca2i = ca2i0
        tau = (2*F*d)/(1e4)
}

BREAKPOINT {
       SOLVE conc METHOD derivimplicit
}

DERIVATIVE conc {
	
	ca2_pump_i = 2*beta_Q10*(ca2i-ca2i0)
	ic = ca2_pump_i
	
	ca2i' = -ica2/tau - beta_Q10*(ca2i-ca2i0)
}
UNITS {
	(mM) = (milli/liter)
      (umho) = (micromho)
        (nA) = (nanoamp)
        (mA) = (milliamp)
        (mV) = (millivolt)
}

NEURON {
	POINT_PROCESS abbott_nmda
	POINTER vpre
	NONSPECIFIC_CURRENT i
	RANGE d, s
	RANGE tauD, tauS, taug, G
	RANGE erev, thresh
}

STATE {
	D 
	S
	g 		(umho)		
	h 		(umho)		
}

PARAMETER {
	v    (mV)
	erev = 0 (mV)
	thresh = 0.5 (mV) 
	d = 0.4
	s = 1.0
	tauD = 300 (ms)
	tauS = 20e3 (ms)
	taug = 2 (ms)
	G    = 10e-12 (umho)
	eta     = 0.33  (/mM)
	mag     = 1     (mM)
	gamma   = 0.06  (/mV)
}	

ASSIGNED {

	i     (nA)
	firing
	vpre  (mV)
}

INITIAL {
	firing = 0
	D = 1
	S = 1
	g = 0
}

BREAKPOINT {

   SOLVE depression METHOD euler 
   aux()
}

DERIVATIVE depression {
   check()
   D' = (1/tauD) * ( 1 - D ) 
   S' = (1/tauS) * ( 1 - S )
   g' = (1/taug) * (-g) 

}

PROCEDURE aux() {
	h = 1/(1 + eta * mag * exp( - (gamma * v)))
	i = (g * h * (v - erev))
}

PROCEDURE check() {




	if (firing && (vpre < thresh)) {
		firing = 0
		}
	if ((vpre >= thresh) && !firing) {
		firing = 1
		D = d * D
		S = s * S
		g = g + G * D * S

		}


}
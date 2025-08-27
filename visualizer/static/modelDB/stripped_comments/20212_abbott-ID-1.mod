UNITS {
      (umho) = (micromho)
        (nA) = (nanoamp)
        (mA) = (milliamp)
        (mV) = (millivolt)
}

NEURON {
	POINT_PROCESS abbott
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
	i = (g * (v - erev))
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

PROCEDURE check2() {
VERBATIM
	if ((v >= thresh) && !firing) {
		firing = 1;
	}
	if (firing && (v < thresh)) {
		firing = 0;
	}
	return 0;
ENDVERBATIM
}
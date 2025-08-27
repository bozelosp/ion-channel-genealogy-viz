UNITS {
	(mV) = (millivolt)
	(mA) = (milliamp)
}

INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}

NEURON {
	SUFFIX InhIAF 
   NONSPECIFIC_CURRENT i
	GLOBAL spikedur, refact, tauAHP, eAHP, gAHPbar
	RANGE Thr, lastspike
	RANGE gPAS, ePAS, gAHP, AHPon, gON, gOFF, eON, eOFF
}


PARAMETER {
	
	gPAS = 0.001		(mho/cm2)
	ePAS = -60			(mV)

	spikedur = 1	(ms)
	refact   = 2.0	(ms)
	Thr	   = -40	(mv)

	gONconst  = 1 	(mho/cm2)
	gOFFconst = 1  (mho/cm2)
	eON  = 40		(mV)
	eOFF = -60		(mV)

	tauAHP   = 0.1 	(/ms)		
	gAHPbar = 0.00005 (mho/cm2)	
	eAHP    = -90 	(mv)
}

ASSIGNED {
	v
	i		(mA/cm2)
	lastspike
	gAHP 		(mho/cm2)
	AHPon				

	gON		(mho/cm2)
	gOFF		(mho/cm2)
}


INITIAL {
	gAHP = 0
	AHPon	  = -9e4

	gON = 0
	gOFF = 0

	lastspike = -9e4
}

BREAKPOINT {
	SOLVE update
	i = gPAS*(v-ePAS) + gAHP*(v-eAHP)+gON*(v-eON)+gOFF*(v-eOFF)
}

PROCEDURE update() { LOCAL q, dv

   gON = 0
   gOFF = 0
   q = (t-lastspike) - spikedur 

   if (q>refact) {				
		if (v>Thr) {				
			gON = gONconst			
			lastspike = t
		}
	}
	else if ( q < 0 ) {			
		gON=gONconst			
	} 
	else if (v > 0) {				
		gOFF = gOFFconst
		gAHP = gAHP + gAHPbar
		AHPon = t
	}
	gAHP = gAHP - gAHP*tauAHP*dt

	VERBATIM
		return 0;
	ENDVERBATIM

}
UNITS {
        (mA) = (milliamp)
        (mV) = (millivolt)
        (S)  = (siemens)
}
 
NEURON {
        SUFFIX nap
        USEION na READ ena WRITE ina
        RANGE  gnabar, ina
}
 
PARAMETER {
	gnabar   =   4e-5 (S/cm2)	

	mvhalf = -52.6		(mV)	
	mslope = -4.6		(mV)	
	mshift = 6.6594 		(mV)

	hvhalf = -48.8		(mV)	
	hslope = 10.0		(mV)	
	hshift = 0.52039			(mV)	

	qfact = 0.36
	
	mtau_base1 = 0.025	
	mtaux1 = 0.14
	mtau_half = -40
	mtau_slope = 10
	mtau_base2 = 0.02
	mtaux2 = 0.145
	
	
}
 
STATE { m h }
 
ASSIGNED {
	ena		(mV)
        v 		(mV)
        ina		(mA/cm2)
        gna		(S/cm2)

        minf
	hinf	

	mtau	(ms)			
   }
 
BREAKPOINT {
        SOLVE state METHOD cnexp
        gna = gnabar * m * h  
        ina = gna * ( v - ena )

}
 

 
INITIAL {
	rates(v)
	
	m = minf
	h = hinf
}

FUNCTION_TABLE tauhnap (v(mV))  (ms)		

DERIVATIVE state { 
        rates(v)
        m' = (minf - m) / mtau
        h' = (hinf - h) / (tauhnap(v)/qfact)    
}
 
PROCEDURE rates(v (mV)) {  
	TABLE minf, hinf, mtau
		FROM -200 TO 200 WITH 201

		minf = 1 / (1 + exp( (v - mvhalf - mshift) / mslope))
		hinf = 1 / (1 + exp( (v - hvhalf - hshift) / hslope))
		
		UNITSOFF
		if (v < mtau_half) {			
			mtau = mtau_base1 + mtaux1 * exp( (v - mtau_half ) / mtau_slope)
		} else {
			mtau = mtau_base2 + mtaux2 * exp( (-v + mtau_half) / mtau_slope)
		}
		UNITSON
}
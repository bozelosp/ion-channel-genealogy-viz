UNITS {
        (mA) = (milliamp)
        (mV) = (millivolt)
        (S)  = (siemens)
}
 
NEURON {
        SUFFIX nap
        USEION na READ ena WRITE ina
        RANGE  gnabar
}
 
PARAMETER {
	gnabar   =   4e-5 (S/cm2)	

	mvhalf = -52.6		(mV)	
	mslope = -4.6		(mV)	

	hvhalf = -48.8		(mV)	
	hslope = 10.0		(mV)	

	qfact = 3
}
 
STATE { m h }
 
ASSIGNED {
	ena		(mV)
        v 		(mV)
        ina		(mA/cm2)
        gna		(S/cm2)

        minf
	hinf	

	taum	(ms)			
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

FUNCTION_TABLE tauh(v(mV))  (ms)		

DERIVATIVE state { 
        rates(v)
        m' = (minf - m) / taum
        h' = (hinf - h) / (tauh(v)/qfact)    
}
 
PROCEDURE rates(v (mV)) {  
	TABLE minf, hinf, taum
		FROM -200 TO 200 WITH 201

		minf = 1 / (1 + exp( (v - mvhalf) / mslope))
		hinf = 1 / (1 + exp( (v - hvhalf) / hslope))
		
		UNITSOFF
		if (v < -40) {			
			taum = 0.025 + 0.14 * exp( (v + 40 ) / 10)
		} else {
			taum = 0.02 + 0.145 * exp( (-v - 40) / 10)
		}
		UNITSON
}
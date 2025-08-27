UNITS {
        (mV) = (millivolt)
	(mA) = (milliamp)      
	(S)  = (siemens)
}
NEURON {
        SUFFIX NaP
        USEION na READ ena WRITE ina
        RANGE  gmax, ina
}
 
PARAMETER {
	gmax   =   1.3802e-7 (S/cm2)	

	mvhalf = -52.6		(mV)	
	mslope = -4.6		(mV)	

	hvhalf = -48.8		(mV)	
	hslope = 10.0		(mV)	
}
 
STATE { m h }
 
ASSIGNED {
	ena		(mV)
        v 		(mV)
        ina		(mA/cm2)
        g		(S/cm2)

        minf
	hinf	

	mtau	(ms)			
   }
 
BREAKPOINT {
        SOLVE state METHOD cnexp
        g = gmax * m * h  
        ina = g * ( v - ena )

}
 

 
INITIAL {
	rates(v)
	m = minf
	h = hinf
}

FUNCTION_TABLE tabhtau (v(mV)) (ms)		

DERIVATIVE state { 
        rates(v)
        m' = (minf - m) / mtau
        h' = (hinf - h) / tabhtau(v)    
}
 
PROCEDURE rates(v (mV)) {  
	TABLE minf, hinf, mtau
		FROM -120 TO 40 WITH 160

		minf = 1 / (1 + exp( (v - mvhalf) / mslope))
		hinf = 1 / (1 + exp( (v - hvhalf) / hslope))
		
		UNITSOFF
		if (v < -40) {			
			mtau = 0.025 + 0.14 * exp( (v + 40 ) / 10)
		} else {
			mtau = 0.02 + 0.145 * exp( (-v - 40) / 10)
		}
		UNITSON
}
UNITS {
        (mA) = (milliamp)
        (mV) = (millivolt)
        (S)  = (siemens)
}
 
NEURON {
        SUFFIX naf
        USEION na READ ena WRITE ina
        RANGE  gnabar, ina, mshift, hshift
		POINTER mu
}
 
PARAMETER {
    gnabar   =   1.5 	(S/cm2)	

	mvhalf = -23.9		(mV)	
	mslope = -11.8		(mV)	
	mshift = 0		(mV)	

	hvhalf = -62.9		(mV)	
	hslope = 10.7		(mV)	
	hshift = 0		(mV)	

	mqfact = 3
	hqfact = 3	
}
 
STATE { m h }
 
ASSIGNED {
		ena				(mV)
        v 				(mV)
        ina 				(mA/cm2)
        gna				(S/cm2)
        minf 
		hinf
		mu (1)
}
 
BREAKPOINT {
        SOLVE state METHOD cnexp
        gna = gnabar * m * m * m  * h
        ina = gna * ( v - ena ) * (1-(mu-1)*0.05)
}
 
 
INITIAL {
	rates(v)
	
	m = minf
	h = hinf
}

FUNCTION_TABLE taumnaf (v(mV))  (ms)	
FUNCTION_TABLE tauhnaf (v(mV))  (ms)	

DERIVATIVE state { 
        rates(v)
        m' = (minf - m) / (taumnaf(v)/mqfact)
        h' = (hinf - h) / (tauhnaf(v)/hqfact)
}
 
PROCEDURE rates(v (mV)) {  
	TABLE minf, hinf DEPEND mshift, hshift, mslope, hslope
		FROM -200 TO 200 WITH 201
			minf = 1 / (1 + exp( (v-mvhalf-mshift) / mslope ) ) 
		    hinf = 1 / (1 + exp( (v-hvhalf-hshift) / hslope ) )
}
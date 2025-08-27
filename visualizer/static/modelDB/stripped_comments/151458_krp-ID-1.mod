UNITS {
        (mA) = (milliamp)
        (mV) = (millivolt)
        (S)  = (siemens)
}
 
NEURON {
        SUFFIX krp
        USEION k READ ek WRITE ik
        RANGE  gkbar, ik
}
 
PARAMETER {
	gkbar   =   0.002 (S/cm2)

	mvhalf = -13.5		(mV)	
	mslope = -11.8		(mV)	
	mshift = 0		(mV)

	hvhalf = -54.7		(mV)	
	hslope = 18.6		(mV)	
 	hshift = 0		(mV)

 	a = 0.7				
 	qfact = 3.0
}
 
STATE { m h }
 
ASSIGNED {
	ek				(mV)
        v 				(mV)
        ik 				(mA/cm2)
        gk				(S/cm2)
        minf 
	hinf
    }
 
BREAKPOINT {
        SOLVE state METHOD cnexp
        gk = gkbar * m * (a*h + (1-a)) 
        ik = gk * ( v - ek )
}
 

 
INITIAL {
	rates(v)
	
	m = minf
	h = hinf
}

FUNCTION_TABLE taumkrp (v(mV))  (ms)
FUNCTION_TABLE tauhkrp (v(mV))  (ms)

DERIVATIVE state { 
        rates(v)
        m' = (minf - m) / (taumkrp(v)/qfact)
        h' = (hinf - h) / (tauhkrp(v)/qfact)
}
 
PROCEDURE rates(v (mV)) {  
	TABLE minf, hinf DEPEND mshift, hshift
		FROM -200 TO 200 WITH 201
			minf = 1 / (1 + exp( (v - mvhalf - mshift) / mslope ))
			hinf = 1 / (1 + exp( (v - hvhalf - hshift) / hslope ))
}
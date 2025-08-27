UNITS {
        (mA) = (milliamp)
        (mV) = (millivolt)
        (S)  = (siemens)
}
 
NEURON {
        SUFFIX kas
        USEION k READ ek WRITE ik
        RANGE  gkbar, ik
}
 
PARAMETER {
    gkbar   =   0.01 (mho/cm2)	

	qfact = 9				
							

	vmh = -27.0	(mV)		
	vmc = -16	(mV)		
	
	vhh = -33.5	(mV)		
	vhc = 21.5	(mV)		
	
	taum0 = 3.4	(ms)		
	Cm = 89.2	(ms)		
	vthm = -34.3	(mV)	
	vtcm = 30.1	(mV)		
	
	alpha = 1				
	vth1 = -0.96	(mV)	
	vtc1 = 29.01	(mV)	
	
	beta = 1				
	vth2 = -0.96	(mV)	
	vtc2 = 100 	(mV)		
	
	Ch = 9876.6	(ms)		
	a = 0.996				
	hshift = 0		(mV)
	htaushift = -90	(mV)	
}
 
STATE { m h }
 
ASSIGNED {
		ek				(mV)
        v 				(mV)
        ik 				(mA/cm2)
        gk				(S/cm2)
        minf
	hinf
        mtau		(ms)
        htau		(ms)
   }
  
INITIAL {
	settables(v)
	m = minf
	h = hinf

}

BREAKPOINT {
        SOLVE state METHOD cnexp
        gk = gkbar * m * m * (a*h + (1-a)) 
        ik = gk * ( v - ek )
}

DERIVATIVE state { 
        settables(v)
	mtau = mtau / qfact
	htau = htau / qfact
        m' = (minf - m)/mtau
        h' = (hinf - h)/htau
}

PROCEDURE settables( v (mV) ) {
	LOCAL left, right

	TABLE minf, hinf, mtau, htau DEPEND hshift, Ch
		FROM -200 TO 200 WITH 201
		
	  	minf = 1 / (1+(exp( (v - vmh) / vmc )))
  		hinf = 1 / (1+(exp( (v - vhh - hshift) / vhc )))
  		
 		mtau = taum0  +  Cm * exp( - ((v-vthm)/vtcm)^2 )

		left = alpha * exp( -(v-vth1-htaushift)/vtc1 )	
		right = beta * exp( (v-vth2-htaushift)/vtc2 )	
		htau = Ch  /  ( left + right )
}
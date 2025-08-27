UNITS {
        (mA) = (milliamp)
        (mV) = (millivolt)
        (S)  = (siemens)
        (molar) = (1/liter)
        (mM) = (millimolar)
}
 
NEURON {
        SUFFIX kir
        USEION k READ ek WRITE ik
        RANGE  gkbar, ik, mvhalf, mslope, mshift, qfact, mu
	
}
 
PARAMETER {
	gkbar  = 0.00015 		(S/cm2)	

	mvhalf = -52		(mV)	
	mslope = 13		(mV)	
	mshift = 30			(mV)	
						
	qfact = 0.5				
	mu = 1
}
 
STATE { m }
 
ASSIGNED {
		ki				(mM)
		ko				(mM)
        v 				(mV)
        ik 				(mA/cm2)
        gk				(S/cm2)
        minf		
        ek				(mV)
		
	
   }
 
BREAKPOINT {
        SOLVE state METHOD cnexp
        gk = gkbar * m
        ik = ((mu-1)*0.25+1)*gk * ( v - ek )
}
  
INITIAL {
	rates(v)
	m = minf
}

FUNCTION_TABLE taumkir (v(mV))  (ms)		

DERIVATIVE state { 
        rates(v)
        m' = (minf - m) / ( taumkir(v)/qfact )
}
 
PROCEDURE rates( v(mV) ) {  
	TABLE minf DEPEND mvhalf, mshift, mslope
		FROM -200 TO 200 WITH 201
			minf = 1  /  ( 1 + exp( (v - mvhalf + mshift) / mslope) )
}
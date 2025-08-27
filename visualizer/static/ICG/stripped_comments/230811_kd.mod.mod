NEURON {
	SUFFIX Ks
	USEION k READ ek WRITE ik
	RANGE gKsbar, ik, gk
	
}

UNITS {
	(mA) = (milliamp)
	(mV) = (millivolt)
        (mM) = (milli/liter)
	
}
INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}
PARAMETER {
	v (mV)
	dt (ms)
	gKsbar	= 1 (mho/cm2) <0,1e9>
	
}


STATE {
	a b
}


ASSIGNED {
	ik (mA/cm2)
	ainf binf
	atau (ms)
	btau (ms)
	gk (mho/cm2)
	ek  (mV)
	
	
}



INITIAL {
	rate(v)
	a = ainf
	b = binf
}

BREAKPOINT {
	SOLVE states METHOD cnexp
		
	gk = gKsbar * a * b
	
	ik = gk*(v-ek)
	
}

DERIVATIVE states {
	rate(v)
	
	a' = (ainf-a)/atau
	b' = (binf-b)/btau
}
UNITSOFF

PROCEDURE rate(v (mV)) {LOCAL va, vb, vc, vd
	
	
	va = v + 34
	vb = v + 65
	vd = v + 63.6
	

if (fabs(va)<1e-04){ va = va+0.00001 }
	   ainf = 1/(1 + exp(-va/6.5))
	   atau=10 
	   
	   

	

if (fabs(vb)<1e-04){ vb = vb+0.00001 }
	   binf = 1/(1 + exp(vb/6.6))

 
if (fabs(vd)<1e-04){ vd = vd+0.00001 }
	
	btau = 3400
	

}

	
UNITSON
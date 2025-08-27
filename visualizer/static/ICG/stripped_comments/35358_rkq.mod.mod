UNITS { 
	(mV) = (millivolt) 
	(mA) = (milliamp) 
} 

NEURON { 
	SUFFIX rkq
	USEION k READ ek WRITE ik
	USEION ca READ cai
	RANGE gmax, ik, g, i
}

PARAMETER { 
	gmax = 1.0 	(mho/cm2)
	v		(mV) 
	ek 		(mV)  
	cai		(1)
}
 
ASSIGNED { 
  g i
  ik 		(mA/cm2) 
  alpha beta	(/ms)
}
 
STATE {
	m
}

BREAKPOINT { 
	SOLVE states METHOD cnexp
        g = gmax * m
	i = g * ( v - ek ) 
        ik = i
}
 
INITIAL { 
	rates( cai )
	m = alpha / ( alpha + beta )
	m = 0
}
 
DERIVATIVE states { 
	rates( cai )
	m' = alpha * ( 1 - m ) - beta * m 
}

UNITSOFF 

PROCEDURE rates(chi) { 

  if (0.2e-4*cai < 0.01) {
    alpha = 0.2e-4*cai
  } else {
    alpha = 0.01
  }
  beta = 0.001
}

UNITSON
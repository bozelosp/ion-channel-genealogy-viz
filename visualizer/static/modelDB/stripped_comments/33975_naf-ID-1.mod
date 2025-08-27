NEURON { SUFFIX naf }
NEURON {  USEION na WRITE ina }
ASSIGNED { ina }
PARAMETER {
	erev 		= 55.  (mV)
	gmax 		= 0.035    (mho/cm2)
        vrest           = 0.

	exptemp		= 27
	maflag 		= 3
	malphaA 	= -0.1
	malphaB		= -10.
	malphaV0	= -35.
	mbflag 		= 1
	mbetaA 		= 4.
	mbetaB		= -18.
	mbetaV0		= -60.
	mq10		= 5
	mexp 		= 3

	haflag 		= 1
	halphaA 	= 0.07
	halphaB		= -20
	halphaV0	= -58.
	hbflag 		= 2
	hbetaA 		= 1.
	hbetaB		= -10.
	hbetaV0		= -28.
	hq10		= 5
	hexp 		= 1

	cao                (mM)
	cai                (mM)

	celsius			   (degC)
	dt 				   (ms)
	v 			       (mV)

	vmax 		= 100  (mV)
	vmin 		= -100 (mV)
} 






INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}

NEURON {
	RANGE gmax, g, i
	GLOBAL erev, Inf, Tau, vmin, vmax, vrest, qq10
} 

CONSTANT {
	  FARADAY = 96489.0	
	  R= 8.31441		

} 

UNITS {
	(mA) = (milliamp)
	(mV) = (millivolt)
	(umho) = (micromho)
} 

ASSIGNED {
	i (mA/cm^2)		
	g (mho/cm^2)
	Inf[2]		
	Tau[2]		
        qq10[2]
} 

STATE { h }

INITIAL { 
 	mh(v)
	h = Inf[1]
}

BREAKPOINT {

  SOLVE states METHOD cnexp
  mh(v)
  g = gmax * Inf[0]*Inf[0]*Inf[0] * h

  i = g*(v-erev) 
  ina=i
} 









DERIVATIVE states {
	mh(v)
	h' = (-h + Inf[1]) / Tau[1]
 }



PROCEDURE mh (v) {
	LOCAL a, b, j
	TABLE Inf, Tau DEPEND maflag, malphaA, malphaB, malphaV0, mbflag, mbetaA, mbetaB, mbetaV0, exptemp, haflag, halphaA, halphaB, halphaV0, hbflag, hbetaA, hbetaB, hbetaV0, celsius, mq10, hq10, vrest, vmin, vmax  FROM vmin TO vmax WITH 200

	qq10[0] = mq10^((celsius-exptemp)/10.)	
	qq10[1] = hq10^((celsius-exptemp)/10.)	

	
	FROM j = 0 TO 1 {
		a = alpha (v, j)
		b = beta (v, j)

		Inf[j] = a / (a + b)
		Tau[j] = 1. / (a + b) / qq10[j]
		if (hexp==0) { Tau[1] = 1. Inf[1] = 1.}
	}
} 


FUNCTION alpha(v,j) {
  LOCAL flag, A, B, V0
  if (j==1 && hexp==0) {
	  alpha = 0
  } else {

     if (j == 1) {
	  A = halphaA B = halphaB V0 = halphaV0+vrest flag = haflag
     } else {
	  A = malphaA B = malphaB V0 = malphaV0+vrest flag = maflag
     }

     if (flag == 1) { 
	 alpha = A*exp((v-V0)/B)	
     } else if (flag == 2) { 
	 alpha = A/(exp((v-V0)/B)+1)
     } else if (flag == 3) { 
	 if(v == V0) {
           alpha = A*B
         } else {
           alpha = A*(v-V0)/(exp((v-V0)/B)-1) }
     }
}
} 


FUNCTION beta (v,j) {
  LOCAL flag, A, B, V0
  if (j==1 && hexp==0) {
	  beta = 1
  } else {

     if (j == 1) {
	  A = hbetaA B = hbetaB V0 = hbetaV0+vrest flag = hbflag
     } else {
	  A = mbetaA B = mbetaB V0 = mbetaV0+vrest flag = mbflag
     }

    if (flag == 1) { 
	 beta = A*exp((v-V0)/B)
     } else if (flag == 2) { 
	 beta = A/(exp((v-V0)/B)+1)
     } else if (flag == 3) { 
	 if(v == V0) {
            beta = A*B 
         } else {
            beta = A*(v-V0)/(exp((v-V0)/B)-1) }
     }
}
}
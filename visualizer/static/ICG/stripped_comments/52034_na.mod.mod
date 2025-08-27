PARAMETER {
	
	gmax 		= 0.1     (umho)




	mvalence 	= 4.3
	mgamma 		=  0.7
	mbaserate 	=  4.2
	mvhalf 		=  -38
	mbasetau 	=  0.05
	mtemp 		=  37
        mq10            =  3
	mexp 		=  3

	hvalence 	= -6
	hgamma		=  0.5
	hbaserate 	=  0.2
	hvhalf 		=  -42
	hbasetau 	=  0.5
	htemp 		=  37
        hq10            =  3
	hexp 		=  1

	cao                (mM)
	cai                (mM)

	celsius		= 37	     (degC)
	dt 				     (ms)
	v 			         (mV)

	vmax 		= 50     (mV)
	vmin 		= -100   (mV)
} 


PROCEDURE iassign() { i = g * (v - ena) ina=i }





INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}

NEURON {
	SUFFIX na
	USEION na READ ena WRITE ina
	
	
	RANGE gmax, g, i, mbaserate
	GLOBAL Inf, Tau, Mult, Add, vmin, vmax
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
        ena (mV)
	i (mA/cm^2)		
	
	ina (mA/cm^2)		
	
	g (mho/cm^2)
	Inf[2]		
	Tau[2]		
	Mult[2]		
	Add[2]		
} 

STATE { m h }

INITIAL { 
 	mh(v)
	if (usetable==0) {
 	  m = Inf[0] h = Inf[1]
	} else {
 	  m = Add[0]/(1-Mult[0]) h = Add[1]/(1-Mult[1]) 
	}
}

BREAKPOINT {

	LOCAL hexp_val, index, mexp_val

	SOLVE states

	hexp_val = 1
	mexp_val = 1

	
	if (hexp > 0) {
		FROM index=1 TO hexp {
			hexp_val = h * hexp_val
		}
	}

	
	if (mexp > 0) {
		FROM index = 1 TO mexp {
			mexp_val = m * mexp_val
		}
	}

	
	
	g = gmax * mexp_val * hexp_val
	iassign()
} 




















PROCEDURE states() {

	

	mh (v*1(/mV))
	m = m * Mult[0] + Add[0]
	h = h * Mult[1] + Add[1]

	VERBATIM
	return 0;
	ENDVERBATIM	
}



PROCEDURE mh (v) {
	LOCAL a, b, j, mqq10, hqq10
	TABLE Add, Mult DEPEND dt, hbaserate, hbasetau, hexp, hgamma, htemp, hvalence, hvhalf, mbaserate, mbasetau, mexp, mgamma, mtemp, mvalence, mvhalf, celsius, mq10, hq10, vmin, vmax  FROM vmin TO vmax WITH 200

	mqq10 = mq10^((celsius-mtemp)/10.)	
	hqq10 = hq10^((celsius-htemp)/10.)	

	
	FROM j = 0 TO 1 {
		a = alpha (v, j)
		b = beta (v, j)

		Inf[j] = a / (a + b)

		VERBATIM
		switch (_lj) {
			case 0:
		
				if ((Tau[_lj] = 1 / (_la + _lb)) < mbasetau) {
					Tau[_lj] = mbasetau;
				}
				Tau[_lj] = Tau[_lj] / _lmqq10;
				break;
			case 1:
				if ((Tau[_lj] = 1 / (_la + _lb)) < hbasetau) {
					Tau[_lj] = hbasetau;
				}
				Tau[_lj] = Tau[_lj] / _lhqq10;
				if (hexp==0) {
					Tau[_lj] = 1.; }
				break;
		}

		ENDVERBATIM
		Mult[j] = exp(-dt/Tau[j])
		Add[j]  = Inf[j]*(1. - exp(-dt/Tau[j]))
	}
} 


FUNCTION alpha(v,j) {
	if (j == 1) {
	   if (hexp==0) {
	     alpha = 1
	   } else {
             alpha = hbaserate * exp((v - hvhalf) * hvalence * hgamma * FRT(htemp)) }
	} else {
		alpha = mbaserate * exp((v - mvhalf) * mvalence * mgamma * FRT(mtemp))
	}
} 


FUNCTION beta (v,j) {
	if (j == 1) {
	   if (hexp==0) {
                beta = 1
	   } else {
		beta = hbaserate * exp((-v + hvhalf) * hvalence * (1 - hgamma) * FRT(htemp)) }
	} else {
		beta = mbaserate * exp((-v + mvhalf) * mvalence * (1 - mgamma) * FRT(mtemp))
	}
} 


FUNCTION FRT(temperature) {
	FRT = FARADAY * 0.001 / R / (temperature + 273.15)
} 


 FUNCTION ghkca (v) { 
       LOCAL nu, efun

       nu = v*2*FRT(celsius)
       if(fabs(nu) < 1.e-6) {
               efun = 1.- nu/2.
       } else {
               efun = nu/(exp(nu)-1.) }

       ghkca = -FARADAY*2.e-3*efun*(cao - cai*exp(nu))
 }
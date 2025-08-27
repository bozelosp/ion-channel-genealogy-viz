UNITS {
	(mV) = (millivolt)
	(mA) = (milliamp)
}

INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}

NEURON {
	SUFFIX Ih
	NONSPECIFIC_CURRENT i
	RANGE gk
	GLOBAL eh,activate_Q10,Q10,gmaxQ10,rate_k,gmax_k,temp1,temp2,tempb
}

PARAMETER {
        v (mV)
	dt (ms)
	gk    = 0.001 (mho/cm2)
	
	celsius

	activate_Q10 = 1
	Q10 = 2.0
	gmaxQ10 = 2.0
	temp1 = 25.0 (degC)
	temp2 = 35.0 (degC)
	tempb = 35.5 (degC)
}

STATE {
        f
}

ASSIGNED {
        eh (mV) 
        i (mA/cm2)
	finf  
	ftau (ms)
	rate_k
	gmax_k
}

BREAKPOINT {
	SOLVE integrate METHOD cnexp
	i = (gk*gmax_k)*f*(v-eh)
}

UNITSOFF

INITIAL {
	LOCAL ktemp,ktempb,ktemp1,ktemp2
	if (activate_Q10>0) {
	  ktemp  = celsius+273.0
	  ktempb = tempb+273.0
	  ktemp1 = temp1+273.0
	  ktemp2 = temp2+273.0
	  rate_k = exp( log(Q10)*((1/ktempb)-(1/ktemp))/((1/ktemp1)-(1/ktemp2)) )
	  gmax_k = exp( log(gmaxQ10)*((1/ktempb)-(1/ktemp))/((1/ktemp1)-(1/ktemp2)) )
	}else{
	  rate_k = 1.0
	  gmax_k = 1.0
	}
        setinf(v)
	f = finf
}

DERIVATIVE integrate {
        setinf(v)
	f' = (finf - f)/ftau
}

PROCEDURE setinf(v) {
                    
	TABLE finf, ftau DEPEND celsius FROM -100 TO 100 WITH 400

	finf = 1.0/(1+exp((v + 80 )/ 5.5))
        ftau = (1.0/(exp(-15.02 - 0.086*v)+exp(-1.5195 + 0.0701*v))) /rate_k
}

UNITSON
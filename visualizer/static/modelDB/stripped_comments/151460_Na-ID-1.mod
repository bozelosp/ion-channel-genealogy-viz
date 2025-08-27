UNITS {
	(mV) = (millivolt)
	(mA) = (milliamp)
}

INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}

NEURON {
	SUFFIX Na
	USEION na READ nai,ena WRITE ina
	RANGE gna, m, h
	GLOBAL rest,activate_Q10,Q10,gmaxQ10,rate_k,gmax_k,temp1,temp2,tempb
}

PARAMETER {
        v (mV)
	dt (ms)
	gna   = 1.0e-7 (mho/cm2)
	rest  = -60.0 (mV) 
	ena
	nai
	celsius
	
	activate_Q10 = 1
	Q10 = 1.980105147e+00
	gmaxQ10 = 1.980105147e+00
	temp1 = 19.0 (degC)
	temp2 = 29.0 (degC)
	tempb = 23.0 (degC)
}

STATE {
        m h  
}

ASSIGNED { 
        ina (mA/cm2)
	alpham (/ms)
	betam (/ms)
	alphah (/ms)
	betah (/ms)
	rate_k
	gmax_k
}

BREAKPOINT {
	SOLVE states METHOD cnexp
	ina  = (gna*gmax_k)*m*m*h*(v-ena)
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
	  
          
          
	  rate_k = 1.60
	  gmax_k = 1.60
	}
        settables(v)
        m = alpham/(alpham+betam)
        h = alphah/(alphah+betah)

}

DERIVATIVE states {
	settables(v)      
	m' = alpham * (1-m) - betam * m
	h' = alphah * (1-h) - betah * h
}

PROCEDURE settables(v) {  
                      
        LOCAL vadj
        TABLE alpham, betam, alphah, betah DEPEND rest,celsius FROM -100 TO 100 WITH 400
	vadj  = v - rest

		
	alpham = rate_k * 0.2 * vtrap((13.1-vadj),4.0)
        betam =  rate_k * 0.175 * vtrap((vadj-40.1),1.0)

                
        alphah = rate_k * 0.08 * exp((17.0-vadj)/18.0)
        betah = rate_k * 2.5 / (exp((40.0-vadj)/5.0) + 1)
}

FUNCTION vtrap(x,y) {  
        if (fabs(x/y) < 1e-6) {
                vtrap = y*(1 - x/y/2)
        }else{
                vtrap = x/(exp(x/y) - 1)
        }
}
 
UNITSON
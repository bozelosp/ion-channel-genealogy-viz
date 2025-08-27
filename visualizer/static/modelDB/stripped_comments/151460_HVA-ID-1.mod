UNITS {
	(mM) = (milli/liter)
	(mV) = (millivolt)
	(mA) = (milliamp)
     FARADAY = (faraday) (coulomb)
           R = (k-mole) (joule/degC)
}

INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}

NEURON {
	SUFFIX HVA
	USEION ca READ cai,cao,eca WRITE ica
	RANGE gcaN, gcaL, iNCa, iLCa
	GLOBAL inactLtau,inactLmax,activate_Q10,Q10,gmaxQ10,rate_k,gmax_k,temp1,temp2,tempb
}

PARAMETER {
        v (mV)
	dt (ms)
        gcaL  = 0.002 (mho/cm2)
	gcaN  = 0.012 (mho/cm2)
	iNCa  = 0.0 (mA/cm2)
	iLCa  = 0.0 (mA/cm2)
	inactLtau = 1220.0 (ms)
	inactLmax = 5.291291201e-01
	eca
	cai
	cao
	celsius

	activate_Q10 = 1
	Q10 = 1.948259241e+00
	gmaxQ10 = 1.948259241e+00
	temp1 = 20.0 (degC)
	temp2 = 30.0 (degC)
	tempb = 22.0 (degC)
}

STATE {
        q u h
}

ASSIGNED { 
        ica (mA/cm2)
	qinf
	uinf
	hinf
	qtau (ms)
	utau (ms)
	htau (ms)
	rate_k
	gmax_k
}

BREAKPOINT {
	LOCAL vghk
	SOLVE states METHOD cnexp
	vghk = ghkg(v,cai,cao,2)
	iNCa = gmax_k*(gcaN * u)*q*q*vghk
	iLCa = gmax_k*(gcaL)*q*q*h*vghk
	ica  = iNCa + iLCa
}


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

        settables(v)
	q = qinf
	u = uinf
	setCadepLinact(cai)
	h = hinf
}

DERIVATIVE states {  
	settables(v)  
	q' = (qinf-q)/qtau
	u' = (uinf-u)/utau
	setCadepLinact(cai)
	h' = (hinf-h)/htau
}

PROCEDURE settables(v) {  
                          
			  
        TABLE qinf, qtau, uinf, utau DEPEND celsius FROM -100 TO 100 WITH 400

                
        qinf   = 1.0/(1.0 + exp((-16.3547869 - v)/11.3))
        qtau   = (1.25/(cosh(-0.031 * (v + 28.8547869)))) /rate_k

                
        uinf   = 1.0/(1.0 + exp((v + 45.3326653)/12.5))
        utau   = (98.0 + cosh(0.021*(24.7673347-v))) /rate_k
}

PROCEDURE setCadepLinact(cai) { 
                
	hinf   = inactLmax+((1.0-inactLmax)/(1.0 + exp((cai-0.7)/0.15)))
        htau   = inactLtau /rate_k
}

INCLUDE "ghk.inc"
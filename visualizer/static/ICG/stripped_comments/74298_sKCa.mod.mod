NEURON {
	SUFFIX sKCa
	USEION ca READ cai
	USEION k READ ki,ek WRITE ik
	RANGE  gk,isKCa
	GLOBAL sKCatau,activate_Q10,Q10,gmaxQ10,rate_k,gmax_k,temp1,temp2,tempb
}

UNITS {
	(mM) = (milli/liter)
	(mA) = (milliamp)
	F = (faraday) (coulombs)	
}

PARAMETER {
        v (mV)
	dt (ms)
	gk = 0.0001 (mho/cm2)
        
	sKCatau = 2.365325544e+01 (ms)
	
	ki
	
	celsius
	
	activate_Q10 = 1
	Q10 = 1.5
	gmaxQ10 = 1.5
	temp1 = 19.0 (degC)
	temp2 = 29.0 (degC)
	tempb = 23.0 (degC)
}

ASSIGNED {
        ek (mV)
        cai (mM)
	ica (mA/cm2)
        ik (mA/cm2)
        winf 
	wtau (ms)
	rate_k
	gmax_k
}

STATE {
        w
}

BREAKPOINT {
	SOLVE integrate METHOD cnexp
	ik = (gk*gmax_k)*w*(v-ek)
	
}


UNITSOFF

INITIAL {
	LOCAL ktemp,ktempb,ktemp1,ktemp2
	if (activate_Q10>0) {
	  rate_k = Q10^((celsius-tempb)/10)
          gmax_k = gmaxQ10^((celsius-tempb)/10)
	}else{
	  rate_k = 1.0
	  gmax_k = 1.0
	}
	setinf(cai)
	w = winf
}

DERIVATIVE integrate {
        setinf(cai)
	w' = (winf - w)/wtau
}

PROCEDURE setinf(cai) {
	LOCAL wcai
	
	wcai = cai*1000
	winf = 0.81/(1+exp((llog(wcai)+0.3)/ -0.46))
	wtau = sKCatau/rate_k
}

FUNCTION llog(x) {  
        if (x>1e-11) {
                llog = log(x)
	}else{
	        llog=0
        }
}

UNITSON
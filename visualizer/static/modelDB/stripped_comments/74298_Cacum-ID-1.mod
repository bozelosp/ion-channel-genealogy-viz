NEURON {
	SUFFIX Cacum
	USEION ca READ ica WRITE cai
	GLOBAL con,cai0,buftau,activate_Q10,Q10,rate_k,temp1,temp2,tempb,depth
}

UNITS {
	(mM) = (milli/liter)
	(mA) = (milliamp)
	F = (faraday) (coulombs)	
}

PARAMETER {
        v (mV)
	dt (ms)
	con   = 0.0			
        Avo   = 6.02e23			
	elc   = 1.602e-19 (coulombs)	
	depth = 200.0 (nm)		
	cai0  = 0.0001(mM)		
	buftau = 1.857456645e+02 (ms)
	cai0_ca_ion
	celsius

	activate_Q10 = 1
	Q10 = 1.2
	temp1 = 19.0 (degC)
	temp2 = 29.0 (degC)
	tempb = 23.0 (degC)
}

ASSIGNED {
	ica (mA/cm2)
        tau (ms)
	rate_k
}

STATE {
	cai (mM)
}

BREAKPOINT {
	SOLVE integrate METHOD cnexp
}

UNITSOFF

INITIAL {
	LOCAL ktemp,ktempb,ktemp1,ktemp2
	if (activate_Q10>0) {
	  rate_k = Q10^((celsius-tempb)/10)
	}else{
	  rate_k = 1.0
	}

	con=1e7/(depth*2.0*Avo*elc)

	tau=buftau/rate_k
	cai=cai0
}

DERIVATIVE integrate {
	cai' = -ica*con + (cai0 - cai)/tau
}

UNITSON
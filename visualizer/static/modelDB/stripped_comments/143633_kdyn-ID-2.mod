NEURON {
	SUFFIX kdyn
	USEION k READ ik WRITE ko, ki 
	RANGE ko, ki, tk, tk0, dep, tkfac
}

UNITS {
	(mM) = (milli/liter)
	(mA) = (milliamp)
	F    = (faraday) (coul)
}

PARAMETER {
	tk0    = 1342.4  (ms)           
	tkfac = 0.025  
	koinf = 3 	(mM)      
	dep   = 290e-3 (micron)     
	KAF   = 1 ()		  
}

ASSIGNED {
	ik     (mA/cm2)
	tk	   (ms)
	ki     (mM)	  
}

STATE { ko (mM) }

BREAKPOINT { 
	SOLVE states METHOD cnexp
}

DERIVATIVE states {
    evaluate_tau()
        ko'= (1e4)*KAF*ik/(F*dep) + (koinf-ko)/tk
}

INITIAL {
	tk = tk0
	ko = koinf
	ki = 106
}

UNITSOFF
PROCEDURE evaluate_tau() { LOCAL tauk 
	tauk = (tkfac*((ko-koinf)^2)) + (tk0/1000) 	
												
	tk = tauk*1000 
	ki = 106 
}
UNITSON
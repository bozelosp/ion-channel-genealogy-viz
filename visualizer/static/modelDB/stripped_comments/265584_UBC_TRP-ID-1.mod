NEURON { 
	SUFFIX Ubc_TRP 
	NONSPECIFIC_CURRENT itrp
	
	USEION ca READ cai
	USEION con_2c READ con_2ci VALENCE 1
	RANGE etrp, gtrp, i, itrp, fcAMP, TonicTRP
} 
 
UNITS { 
	(mA) = (milliamp) 
	(mV) = (millivolt) 
} 
 
PARAMETER { 
	v (mV) 
	gtrp = 4.18e-6 (mho/cm2)
	celsius = 30 (degC)
	etrp = 0 (mV)
	fcAMP = 1
	theta = 0 (1)
	TonicTRP = 0.05 
    } 
    
    ASSIGNED { 
	cAMPi  (mM)
	cai (mM)
	itrp (mA/cm2) 
	i (mA/cm2)
	con_2ci (mM)
    }
    
    BREAKPOINT { 
	
	if(con_2ci > 0.0075){
	     itrp = (gtrp*(v - etrp) * (1 + (3*con_2ci)))
         }else {
          itrp = 0    
	}
	i = itrp
	
}
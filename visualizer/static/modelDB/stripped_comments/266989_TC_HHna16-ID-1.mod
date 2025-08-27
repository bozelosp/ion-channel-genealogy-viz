INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}

NEURON {
	SUFFIX TC_HHna16
	USEION na READ ena WRITE ina
	
	
        RANGE gna_max, vtraub, i_rec
	RANGE m_inf, h_inf 
	RANGE tau_m, tau_h 
	
	RANGE ina 
}


UNITS {
	(mA) = (milliamp)
	(mV) = (millivolt)
	(S)  = (siemens)
}

PARAMETER {
	gna_max	= 1.0e-1 	(S/cm2) 


	celsius         (degC)
	dt              (ms)
	v               (mV)
	vtraub = -55.5  

}

STATE {

  m h 
}

ASSIGNED {
	ina	(mA/cm2)
	
	ena	(mV)
	
	i_rec	(mA/cm2)
	m_inf
	h_inf
	
	tau_m
	tau_h
	
	
	
	
	tcorr
}


BREAKPOINT {
	SOLVE states METHOD cnexp
	ina   = gna_max * m*m*m*h * (v - ena)
	
	
        i_rec = ina
}


DERIVATIVE states {   
	evaluate_fct(v)
	m' = (m_inf - m) / tau_m
	h' = (h_inf - h) / tau_h
	
}









UNITSOFF
INITIAL {





  	tcorr = 3.0 ^ ((celsius-36)/ 10 )
        evaluate_fct(v)
	m = m_inf
	h = h_inf


}

PROCEDURE evaluate_fct(v(mV)) { LOCAL a,b,v2, v3

	v2 = v - vtraub + 15
        v3 = v - vtraub + 15
	

	if(v2 == 13 || v2 == 40 || v2 == 15 ){
    	v = v+0.0001
    }

	a = 0.32 * (13-v2) / ( exp((13-v2)/(1.2*4)) - 1)
	b = 0.28 * (v2-40) / ( exp((v2-40)/(1.2*5)) - 1)
	tau_m = 1 / (a + b) / tcorr
	m_inf = a / (a + b)

                                
        a = 0.128 * exp((17-v3)/18) 
        b = 4 / ( 1 + exp((40-v3)/5) )
	tau_h = 1 / (a + b) / tcorr
        h_inf = a / (a + b)

                                
 

	
	
	
	

	
	
	

}

UNITSON
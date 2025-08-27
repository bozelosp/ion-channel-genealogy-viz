INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}

NEURON {
        SUFFIX Golgi_Ca_LVA
        
        USEION ca READ eca,cai,cao WRITE ica
        RANGE g, gca2bar, m_inf, tau_m, h_inf, tau_h, shift
	RANGE m ,h
	RANGE phi_m, phi_h
	RANGE v0_m_inf,v0_h_inf,k_m_inf,k_h_inf,C_tau_m
	RANGE A_tau_m,v0_tau_m1,v0_tau_m2,k_tau_m1,k_tau_m2
	RANGE C_tau_h ,A_tau_h ,v0_tau_h1,v0_tau_h2,k_tau_h1 ,k_tau_h2

    }

UNITS {
        (molar) = (1/liter)
        (mV) =  (millivolt)
        (mA) =  (milliamp)
        (mM) =  (millimolar)

        FARADAY = (faraday) (coulomb)
        R = (k-mole) (joule/degC)
}

PARAMETER {
        v               (mV)
        celsius (degC)
        
	   gca2bar  = 2.5e-4 (mho/cm2)
        shift   = 2     (mV)            
        cai  (mM)           
        cao  (mM)
	
	v0_m_inf = -50 (mV)
	v0_h_inf = -78 (mV)
	k_m_inf = -7.4 (mV)
	k_h_inf = 5.0  (mv)
	
	C_tau_m = 3
	A_tau_m = 1.0
	v0_tau_m1 = -25 (mV)
	v0_tau_m2 = -100 (mV)
	k_tau_m1 = 10 (mV)
	k_tau_m2 = -15 (mV)
	
	C_tau_h = 85
	A_tau_h = 1.0
	v0_tau_h1 = -46 (mV)
	v0_tau_h2 = -405 (mV)
	k_tau_h1 = 4 (mV)
	k_tau_h2 = -50 (mV)
	
    }
    

STATE {
        m h
}

ASSIGNED {
        ica     (mA/cm2)
        eca     (mV)
	
	g        (mho/cm2) 
        m_inf
        tau_m   (ms)
        h_inf
        tau_h   (ms)
        phi_m
        phi_h
}

BREAKPOINT {
        SOLVE ca2state METHOD cnexp
        
        g = gca2bar * m*m*h
        ica = gca2bar * m*m*h * (v-eca)
}

DERIVATIVE ca2state {
        evaluate_fct(v)

        m' = (m_inf - m) / tau_m
        h' = (h_inf - h) / tau_h
}

UNITSOFF
INITIAL {







        evaluate_fct(v)
        m = m_inf
        h = h_inf
}

PROCEDURE evaluate_fct(v(mV)) { 



        phi_m = 5.0 ^ ((celsius-24)/10)
        phi_h = 3.0 ^ ((celsius-24)/10)
	
	TABLE m_inf, tau_m, h_inf, tau_h
	DEPEND shift, phi_m, phi_h FROM -100 TO 30 WITH 13000 
        m_inf = 1.0 / ( 1 + exp((v + shift - v0_m_inf)/k_m_inf) )
        h_inf = 1.0 / ( 1 + exp((v + shift - v0_h_inf)/k_h_inf) )
	
        tau_m = ( C_tau_m + A_tau_m / ( exp((v+shift - v0_tau_m1)/ k_tau_m1) + exp((v+shift - v0_tau_m2)/k_tau_m2) ) ) / phi_m
        tau_h = ( C_tau_h + A_tau_h / ( exp((v+shift - v0_tau_h1)/k_tau_h1) + exp((v+shift - v0_tau_h2)/k_tau_h2) ) ) / phi_h
}
UNITSON
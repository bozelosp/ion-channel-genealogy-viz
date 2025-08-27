INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}

NEURON {
	SUFFIX ia
	USEION k  READ ek WRITE ik VALENCE 1
	RANGE gmax, i
	RANGE m_inf, tau_m, h_inf, tau_h, exptemp, q10
}

UNITS {
	(mV) =	(millivolt)
	(mA) =	(milliamp)
}

PARAMETER {
  ek
  v		(mV)
  
  gmax	= 0.0	(mho/cm2)
  exptemp= 23.5
  q10 = 3
}

STATE {
  m h
}

ASSIGNED {
	ik	(mA/cm2)
	i	(mA/cm2)
	m_inf
	tau_m	(ms)
	h_inf
	tau_h	(ms)
        tadj
}

BREAKPOINT {
	SOLVE states METHOD cnexp
	i = gmax * (m*m*m*m*h * (v-ek))
        ik = i
}

DERIVATIVE states {
	evaluate_fct(v)

	m' = (m_inf - m) / tau_m
	h' = (h_inf - h) / tau_h 
}

UNITSOFF
INITIAL {
        tadj = pow(q10,((celsius-exptemp)/10))
	evaluate_fct(v)
	m = m_inf
	h = h_inf




}

PROCEDURE evaluate_fct(v(mV)) { 
  
  

  m_inf = 1.0 / ( 1 + exp(-(v+60)/8.5) )
  h_inf = 1.0 / ( 1 + exp((v+78)/6.0) )

  tau_m = (1.0/  (exp((v+35.82)/19.69)+exp(-(v+79.69)/12.7)) +0.37) / tadj

  if (v < -63) {
    tau_h =  1.0 /(exp((v+46.05)/5)+exp(-(v+238.4)/37.45)) / tadj

  } else {	
    tau_h = 19.0/tadj
   
  }
}
UNITSON
NEURON {
       SUFFIX akca
       USEION k READ ek WRITE ik
	   USEION ca READ cai
       RANGE gbar, ek, ik

}

UNITS {
      (S) = (siemens)
      (mV) = (millivolts)
      (mA) = (milliamp)
	  (molar) = (/liter)
	  (mM) = (millimolar)
}

PARAMETER {
    gbar =0.00022989 (S/cm2)
	Q10kcac=2.30

	
	A_alphac=750.0 (/ms-mM) 
	B_alphac=-10.0 (mV)
	C_alphac=12.0	(mV)
	
	A_betac=0.05 (/ms)	
	B_betac=-10.0 (mV)
	C_betac=-60.0 (mV)
}


ASSIGNED {
	 v	(mV) 
	 ik	(mA/cm2)
	 cai (mM)
	 celsius (degC)
	 g	(S/cm2)
	 tau_c	(ms)
	 cinf
	 alphac (/ms)
	 betac (/ms)
     ek	(mV)
         
}

STATE { c }

BREAKPOINT {
	   SOLVE states METHOD cnexp
	   g = gbar * c
	   ik = g * (v-ek)
}

INITIAL {
	rates(v) 
	
	
    c = cinf
}

DERIVATIVE states {
	   rates(v)
	   c' = (cinf - c)/tau_c
}

FUNCTION alpha(Vm (mV)) (/ms) {
	alphac=A_alphac*cai*exp((Vm+B_alphac)/C_alphac)
}

FUNCTION beta(Vm (mV)) (/ms) {
	betac=A_betac*exp((Vm+B_betac)/C_betac)
}

FUNCTION rates(Vm (mV)) (/ms) {
	alpha(Vm)
	beta(Vm)
	
	tau_c = 4.5/(alphac+betac)
	cinf = alphac/(alphac+betac)
	
	if (celsius >= 37) {
	tau_c=tau_c*(1/Q10kcac)
	}
}
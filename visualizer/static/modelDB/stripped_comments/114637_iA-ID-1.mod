INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}

NEURON {
	SUFFIX iA
	USEION k READ ek WRITE ik
        RANGE gkbar, m_inf1, tau_m, h_inf, tau_h1, ik, actvha, inactvha
}

UNITS {
	(mA) = (milliamp)
	(mV) = (millivolt)
}

PARAMETER {
	v		(mV)
	celsius		(degC)
	dt		(ms)
	ek		(mV)
	gkbar= 0.00345	(mho/cm2)
	actvha=50 	(mV)
	inactvha=60	(mV)
	
}

STATE {
	m1 h1
}

ASSIGNED {
	ik		(mA/cm2)
	m_inf1[1]
	tau_m[1]		(ms)
	h_inf[1]
	tau_h1[1]	(ms)
	tadj
}

BREAKPOINT {
	SOLVE states 
 	ik = gkbar * m1^4*h1 * (v-ek)
}








PROCEDURE states() {
        
	mhn(v)
	m1= m1 + (1-exp(-dt/tau_m[0]))*(m_inf1[0]-m1)
	h1= h1 + (1-exp(-dt/tau_h1[0]))*(h_inf[0]-h1)
}


INITIAL {



}




	
	
	







	
	
	




FUNCTION vartauh(v,i) {LOCAL a,b
	if (i==0) {
		if (v<-63) {vartauh =1.0 / ((exp((v+46.05)/5)+exp(-(v+238.4)/37.45))) / tadj
		} else {vartauh = 26.9584/tadj}
	}
}

FUNCTION vartaum(v,i) {LOCAL a,b
	if (i==0) {
		vartaum = (1.0/(exp((v+35.82)/19.69)+exp(-(v+79.69)/12.7))+0.37) / tadj
	}
}

FUNCTION varssm(v, i) {
	if (i==0) {
		varssm = 1.0 / (1+exp(-(v+actvha)/8.5))
	}
}

FUNCTION varssh(v, i) {
	if (i==0) {
		varssh = 1.0/(1+exp((v+inactvha)/6))
	}
}

PROCEDURE mhn(v) {
	TABLE m_inf1,tau_m,h_inf,tau_h1
	DEPEND inactvha,celsius, actvha,dt
	FROM -100 TO 100 WITH 2000 
	FROM i=0 TO 0 {
		m_inf1[i] = varssm(v,i)
		tau_m[i]=vartaum(v,i)
		h_inf[i] = varssh(v,i)
		tau_h1[i]=vartauh(v,i)
	}
	tadj = 3^((celsius-23.5)/10)
}
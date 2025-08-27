INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}

NEURON {
	SUFFIX tcihshift
	NONSPECIFIC_CURRENT ih
	RANGE ghbar, eh
	RANGE q_inf
	RANGE tau_q, vshifto, vsh_increment, tau_vsh, nsm
	POINTER vext, pmodyn	
}

UNITS {
	(mA) = (milliamp)
	(mV) = (millivolt)
}

PARAMETER {
	ghbar = 0.00007 (mho/cm2)
	eh	= -43	(mV)
	qexp	= 3
	vshifto	= 5	(mV)
	celsius		(degC)
	dt              (ms)
	v               (mV)
	vsh_increment = 0.015 (mV/ms)
	tau_vsh = 2000 (ms)
	stimon = 0
	nsm = 1
}

STATE {
	q
	vsh (mV)
}

ASSIGNED {
	ih	(mA/cm2)
	q_inf
	tau_q (ms)
	vext 
	pmodyn
	tadj
}

BREAKPOINT {
	SOLVE states METHOD euler
	ih  = ghbar * q^qexp * (v-eh)
}

DERIVATIVE states {   
       evaluate_fct(v,vsh)
	q' = (q_inf - q) / tau_q
	vsh' = (vshifto - vsh) / tau_vsh - nsm*pmodyn*vsh_increment

}

UNITSOFF
INITIAL {



	tadj = 3.0 ^ ((celsius-36)/ 10 )

	evaluate_fct(v, vsh)
	q = q_inf
	vsh = vshifto
}

PROCEDURE evaluate_fct(v(mV), vsh(mV)) {
	tau_q = 1 / (Exp(-14.59-(0.086*(v+vsh))) + Exp(-1.87+(0.0701*(v+vsh))))
	q_inf = 1 / (Exp((v+vsh+75)/5.5) + 1)
	if(vext==0) {      
		stimon=0
	}else{
		stimon=1
	}
}

FUNCTION Exp(x) {
	if (x < -100) {
		
	}else{
		if (x > 700) {
			Exp = exp(700)
		}else{
			Exp = exp(x)
		}
	}
}
UNITSON
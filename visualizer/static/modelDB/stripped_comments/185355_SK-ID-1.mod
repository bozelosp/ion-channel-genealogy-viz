UNITS {
        (molar) = (1/liter)
        (mM)    = (millimolar)
	(mA)	= (milliamp)
	(mV)	= (millivolt)
}

NEURON {
	SUFFIX sk
	USEION k READ ek WRITE ik 
	USEION nca READ ncai VALENCE 2
	USEION lca READ lcai VALENCE 2
	USEION tca READ tcai VALENCE 2
	RANGE gsk, gskbar, qinf, qtau, ik
}

INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}

PARAMETER {
	celsius=6.3 (degC)
	v		(mV)
	dt		(ms)
	gskbar  (mho/cm2)
	ek	(mV)
	cai (mM)
	ncai (mM)
	lcai (mM)
	tcai (mM)
}

STATE { q }

ASSIGNED {
	ik (mA/cm2) gsk (mho/cm2) qinf qtau (ms) qexp
}


BREAKPOINT {          
	SOLVE state
        gsk = gskbar * q*q
	ik = gsk * (v-ek)
}

UNITSOFF

INITIAL {
	cai = ncai + lcai + tcai	
	rate(cai)
	q=qinf
	VERBATIM
	ncai = _ion_ncai;
	lcai = _ion_lcai;
	tcai = _ion_tcai;
	ENDVERBATIM
}


PROCEDURE state() {  
	cai = ncai + lcai + tcai
	rate(cai)
	q = q + (qinf-q) * qexp
	VERBATIM
	return 0;
	ENDVERBATIM
}

LOCAL q10
PROCEDURE rate(cai) {  
	LOCAL alpha, beta, tinc
	q10 = 3^((celsius - 6.3)/10) 
		

        
	alpha = 0.00246/exp((12*log10(cai)+28.48)/-4.5)
	beta = 0.006/exp((12*log10(cai)+60.4)/35)
	qtau = 1 / (alpha + beta)
	qinf = alpha * qtau
	tinc = -dt*q10
	qexp = 1 - exp(tinc/qtau)*q10
}

UNITSON
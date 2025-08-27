TITLE Hill-Mashma model
 
COMMENT

ENDCOMMENT

NEURON {
    SUFFIX fHill
    RANGE F_norm, Kse, p0_5
	RANGE g1, g2, g3
	RANGE a0, b0, c0, d0
	USEION mg READ mgi VALENCE 2
	USEION cl READ cli
}

PARAMETER {
	Kse = 0.16
	p0_5 = 100.3
	g1 = -0.0045		
	g2 = 0.0981		
	g3 = 0.6211
	a0 = 0.4		
	b0 = 99.7		
	c0 = -57.1		
	d0 = 42.2	
	xm_init = 5		
	xce_init = 5	
}

STATE {
	A
	xce
	xm
}

ASSIGNED {
	F_norm
	Fc_norm
	mgi		
	cli		
}

BREAKPOINT { LOCAL d_xm, d_xce, d_se
	A = mgi		
	xm = cli	

	SOLVE state_hill METHOD cnexp
		
	F_norm = p0_5*Kse*xse(xm, xce)/p0_5
}

DERIVATIVE state_hill {
	Fc_norm = p0_5*g(xm)*A/p0_5
	xce' = dxdt (F_norm, Fc_norm)
}

FUNCTION xse (x, y) { LOCAL d_xm, d_xce, d_se
	d_xm = xm - xm_init
	d_xce = xce - xce_init
	d_se = d_xm - d_xce
	
	if (d_se <= 0) {xse = 0} 
	else {xse = d_se}
}

FUNCTION g (x) {
	:g = exp(-((x-g1)/g2)^2)+g3
	g = g1*x^2+g2*x+g3
}

FUNCTION dxdt (x, xc) {LOCAL gain_length 
	if (x <= xc) {
		dxdt = (10^-3)*(-b0*(xc-x))/(x+a0*xc/p0_5)
	} else {
		gain_length = (-d0*(xc-x))/(2*xc-x+c0*xc/p0_5)
		if (gain_length <= 0) {dxdt = (10^-3)*1e5}
		else {dxdt = (10^-3)*gain_length}
	}
}

INITIAL {
	A = 0
	xm = xm_init
	xce = xce_init
	F_norm=1e-5
}
NEURON {
	SUFFIX hva
	
	USEION ca READ eca WRITE ica
	RANGE gbar
}

UNITS {
	(S)	=	(siemens)
	(mV)	=	(millivolt)
	(mA)	=	(milliamp)
}

PARAMETER {
	gbar = 1110e-6	(S/cm2) < 0, 1e9 >
	
}

ASSIGNED {
        eca (mV)
	ica (mA/cm2)
	
	v (mV)
	g (S/cm2)
	sinf
	rinf
	tau_s (ms)
	tau_r (ms)
}

STATE {	s r }

BREAKPOINT {
	SOLVE states METHOD cnexp
	g = gbar * s*s * r
	ica = g * (v - eca)
	
}

INITIAL {
	
	s = alphas(v)/(alphas(v) + betas(v))
	r = alphar(v)/(alphar(v) + betar(v))
}

DERIVATIVE states {
	rates(v)
	s' = (sinf - s)/tau_s
	r' = (rinf - r)/tau_r
}

LOCAL alpha_r, alpha_s	

FUNCTION alphas(Vm (mV)) (/ms) {
	UNITSOFF
	alphas = 2.0 /(1 + exp(-(Vm + 2.0)* 0.054))
	UNITSON
}

FUNCTION betas(Vm (mV)) (/ms) {
	UNITSOFF
	betas =  -0.08 * (Vm + 15.9) / (1 - exp( (Vm + 15.9)*0.2 ))
	UNITSON
}

FUNCTION taus(Vm (mV)) (/ms) {
	UNITSOFF
	taus = 1.0 / (alpha_s + betas(Vm))	
	UNITSON
}

FUNCTION alphar(Vm (mV)) (/ms) {
	UNITSOFF
	alphar = 0.01 * exp( -(Vm + 60)/20 )
	if (alphar > 0.01) {
		alphar = 0.01
	}
	UNITSON
}

FUNCTION betar(Vm (mV)) (/ms) {
	UNITSOFF
	betar =  0.01 - alphar(Vm)
	UNITSON
}

FUNCTION taur(Vm (mV)) (/ms) {
	UNITSOFF
	taur = 1.0 / (alpha_r + betar(Vm))	
	UNITSON
}





PROCEDURE rates(Vm (mV)) {
	alpha_s = alphas(Vm)
	tau_s = taus(Vm)
	sinf = alphas(Vm) * tau_s

	alpha_r = alphar(Vm)
	tau_r = taur(Vm)
	rinf = alpha_r * tau_r
}
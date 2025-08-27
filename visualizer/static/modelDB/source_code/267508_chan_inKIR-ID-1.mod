TITLE inKIR channel

COMMENT

ENDCOMMENT

NEURON {
	SUFFIX inKIR
	USEION  k READ ek WRITE ik
	RANGE  g, ik
	GLOBAL minf, mtau, hinf, htau, eff_hinf, a
}

UNITS {
	(mA) = (milliamp)
	(uA) = (microamp)
	(mV) = (millivolt)
	(mS) = (millimho)
}

CONSTANT {	
	Q10 = 3 (1)
}

PARAMETER {
	ek			(mV)
	gmax = 1.4e-4	(mho/cm2)	<0,1e9>
	m_vh = -82	(mV)	: half activation
	m_ve = 13		(mV)	: slope
	a = 0.47 (1)		: 0.47 default, 0.27 for pinKir
:	a = .27 (1)		
}

ASSIGNED {
	v	(mV)
	g	(mho/cm2)
	ik	(mA/cm2)
	minf	(1)
	mtau	(ms)
	hinf (1)
	htau (ms)
	eff_hinf (1)
	qt (1)
}

STATE {
	m h
}

BREAKPOINT {
	SOLVE states METHOD cnexp
:	g = gmax*m*h
	g=gmax*m*(a*h+(1 -a)) 
	ik = g*(v - ek)
}

INITIAL {
	qt = Q10^((celsius-35)/10)
	rates(v)
	m = minf
	h = hinf
}

DERIVATIVE states { 
	rates(v)
	m' = (minf-m)/mtau
	h' = (hinf-h)/htau
}

FUNCTION_TABLE tabmtau(v(mV)) (ms)
FUNCTION_TABLE tabhtau(v(mV)) (ms)
FUNCTION_TABLE tabhinf(v(mV))(1)

: rates() computes rate and other constants at present v
: call once from hoc to initialize inf at resting v

PROCEDURE rates(v(mV)) {
	mtau = tabmtau(v)/qt
:	mtau = tabmtau(v)
	htau = tabhtau(v)/qt
:	htau = tabhtau(v)
	hinf = tabhinf(v)
	eff_hinf = a*hinf+(1 -a)
	minf = 1/(1 + exp((v - m_vh)/m_ve))
}


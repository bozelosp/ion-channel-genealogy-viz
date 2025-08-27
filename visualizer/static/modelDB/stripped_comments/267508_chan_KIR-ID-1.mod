NEURON {
	SUFFIX KIR
	USEION  k READ ek WRITE ik
	RANGE g, ik
	GLOBAL minf, mtau, gmax
}

CONSTANT {
	Q10 = 3 (1)
}

UNITS {
	(mA) = (milliamp)
	(uA) = (microamp)
	(mV) = (millivolt)
	(mS) = (millimho)
}

PARAMETER {
	ek			(mV)
	gmax = 1.4e-4	(mho/cm2)	<0,1e9>
	m_vh = -82	(mV)	
	m_ve = 13		(mV)	
}

ASSIGNED {
	v	(mV)
	g	(mho/cm2)
	ik	(mA/cm2)
	minf	(1)
	mtau	(ms)
	qt (1)
}

STATE {
	m
}

BREAKPOINT {
	SOLVE states METHOD cnexp
	g = gmax*m
	ik = g*(v - ek)
}

INITIAL {
	qt = Q10^((celsius-35)/10)
	rates(v)
	m = minf
}

DERIVATIVE states { 
	rates(v)
	m' = (minf-m)/mtau
}

FUNCTION_TABLE tabmtau(v(mV)) (ms)




PROCEDURE rates(v(mV)) {

	mtau = tabmtau(v)/qt
	minf = 1/(1 + exp((v - m_vh)/m_ve))
}
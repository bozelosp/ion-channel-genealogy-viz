NEURON {
	SUFFIX KRP
	USEION  k READ ek WRITE ik
	RANGE g, gmax, ik
	GLOBAL minf, mtau, hinf, htau
}

UNITS {
	(mA) = (milliamp)
	(uA) = (microamp)
	(mV) = (millivolt)
	(mS) = (millimho)
}

PARAMETER {
	gmax = 0.001	(mho/cm2)	<0,1e9>
	ek 		        (mV)
	m_vh = -13.5	(mV)	
	m_ve = -11.8	(mV)	
	h_vh = -54.7	(mV) 
	h_ve = 18.6	(mV) 
}

ASSIGNED {
	v	(mV)
	g	(mho/cm2)
	ik	(mA/cm2)
	minf	(1)
	hinf   (1)
	mtau	(ms)
	htau  (ms)
}

CONSTANT {
	a = 0.7 (1)
}

STATE {
	m
	h
}

BREAKPOINT {
	SOLVE states METHOD cnexp
	g = gmax*m*(a*h+(1-a))
	ik = g*(v - ek)
}

INITIAL {
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



PROCEDURE rates(v(mV)) {

	mtau = tabmtau(v)
	minf = 1/(1 + exp((v - m_vh)/m_ve))
	htau = tabhtau(v)
	hinf = 1/(1 + exp((v - h_vh)/h_ve))
}
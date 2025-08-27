UNITS {
	(molar) = (1/liter)
	(mV) =	(millivolt)
	(mA) =	(milliamp)
	(mM) =	(millimolar)
	FARADAY = (faraday)  (kilocoulombs)
	R = (k-mole) (joule/degC)
	(nA)=(nanoamp)
	(um)=(micrometer)
}

NEURON {
	SUFFIX skca3
	USEION ca READ cai
	USEION k READ ek WRITE ik
	RANGE ik, gbar, g
	RANGE oinf, hcsk3, E50hsk3
	RANGE m_vh, m_sf, m
	THREADSAFE
}

PARAMETER {
	v				(mV)
	gbar = 0.0009	(mho/cm2)
	cai				(mM) 
	ek				(mV)
		
	
	hcsk3	= 5.6 	(1)
	E50hsk3 = 0.42e-3 (mM)

	m				(1)
	m_vh = 24		(mV) 
	m_sf = 128		(mV)
}

ASSIGNED {
	ik		(mA/cm2)
	o
	oinf
	g		(mho/cm2)
}

BREAKPOINT {
	rate(cai,v)
	g = gbar*o*m
	ik = gbar*o*m*(v - ek)
}

INITIAL {
	rate(cai,v)
}

FUNCTION_TABLE tabvh(cai(mM)) (mV) 
FUNCTION_TABLE tabsf(cai(mM)) (mV) 

PROCEDURE rate(ca (mM), v(mV)) {
	UNITSOFF
	
	
	oinf = (ca^hcsk3)/((E50hsk3^hcsk3)+(ca^hcsk3))
	UNITSON
	o=oinf
	
	
	
	m_vh = tabvh(ca)
	m_sf = tabsf(ca)
	
	m = 1/(1+exp((v-(ek+m_vh))/(m_sf)))
}
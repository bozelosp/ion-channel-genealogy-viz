UNITS {
	(mA) =(milliamp)
	(mV) =(millivolt)
	(uF) = (microfarad)
	(molar) = (1/liter)
	(nA) = (nanoamp)
	(mM) = (millimolar)
	(um) = (micron)
	FARADAY = 96520 (coul)
	R = 8.3134	(joule/degC)
}
 
NEURON { 
	SUFFIX ch_Navngf 
	USEION na READ ena WRITE ina VALENCE 1
	RANGE g, gmax, minf, mtau, hinf, htau, ina, m, h
	RANGE myi
	THREADSAFE
}

PARAMETER {
	ena  (mV)
	gmax (mho/cm2)   
	
	mAlphC = -0.34133 (1)
	mAlphV = 24 (mV)
	mBetaC = 0.28483 (1)
	mBetaV = -4 (mV)

	hAlphC = 0.29648 (1)
	hAlphV = 64.4184 (mV)
	hBetaC = 3.0931  (1)
	hBetaV = 12.1463 (mV)
}

STATE {
	m h
}

ASSIGNED {
	v (mV) 
	celsius (degC) 
	dt (ms) 

	g (mho/cm2)
	ina (mA/cm2)
	minf
	hinf
	mtau (ms)
	htau (ms)
	mexp
	hexp 
	myi (mA/cm2)
}

BREAKPOINT {
	SOLVE states
	g = gmax*m*m*m*h  
	ina = g*(v - ena)
	myi = ina
}
 
UNITSOFF
 
INITIAL {
	trates(v)
	m = minf
	h = hinf
}

PROCEDURE states() {	
	trates(v)			
	m = m + mexp*(minf-m)
	h = h + hexp*(hinf-h)
}
 
LOCAL q10	
PROCEDURE rates(v) {  
                      
	LOCAL  alpha, beta, sum	

	q10 = 3^((celsius - 34)/10)

	
	alpha = mAlphC*vtrap((v+mAlphV),-5)
	beta = mBetaC*vtrap((v+mBetaV),5)
	sum = alpha+beta        
	mtau = 1/sum 
	minf = alpha/sum
	
	
	alpha = hAlphC/exp((v+hAlphV)/20)
	beta = hBetaC/(1+exp((v+hBetaV)/-10))
	sum = alpha+beta
	htau = 1/sum 
	hinf = alpha/sum 		
}
 
PROCEDURE trates(v) {  
                      
	LOCAL tinc	
	TABLE minf, mexp, hinf, hexp, mtau, htau
	DEPEND dt, celsius, mAlphV, mAlphC, mBetaV, mBetaC, hAlphV, hAlphC, hBetaV, hBetaC
	FROM -100 TO 100 WITH 200

	rates(v)	
				
				

	tinc = -dt * q10

	mexp = 1 - exp(tinc/mtau)
	hexp = 1 - exp(tinc/htau)
}
 
FUNCTION vtrap(x,y) {  
	if (fabs(x/y) < 1e-6) {
		vtrap = y*(1 - x/y/2)
	}else{  
		vtrap = x/(exp(x/y) - 1)
	}
}
 
UNITSON
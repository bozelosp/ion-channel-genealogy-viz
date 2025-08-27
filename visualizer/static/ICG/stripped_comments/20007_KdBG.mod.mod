NEURON {
	SUFFIX KdBG
	USEION k READ ek WRITE ik
	RANGE  gbar,ik
	GLOBAL xtau, ytau, xinf, yinf
}

UNITS {
	(S)	= (siemens)
	(mA)	= (milliamp)
	(mV)	= (millivolt)
	FARADAY	= (faraday) (coulombs)
	R	= (k-mole)  (joule/degC)
}

PARAMETER {
	gbar	=   1.0e-3	(S/cm2)
	celsius			(degC)
	Ky	=   2.0e-4	(1/ms)
	gammay	=   0.0		(1)
	zettax	=   3.0		(1)
	zettay	=  -2.5		(1)
	vhalfx	= -63.0		(mV)
	vhalfy	= -73.0		(mV)
	taox	=   1.0		(ms)
	taoy	=   0.0		(ms)
}

ASSIGNED {
        ek      (mV)
	v       (mV)
	ik     	(mA/cm2)
	xtau    (ms)
	ytau    (ms)
	xinf	(1)
	yinf	(1)
	q10	(1)
	T     	(K)
}

STATE { xs ys }

BREAKPOINT { 
	SOLVE states METHOD cnexp
	ik= gbar * xs^4 * ys^4 * ( v - ek ) 
}

DERIVATIVE states {
	rates()
	xs'= (xinf- xs)/ xtau	
	ys'= (yinf- ys)/ ytau
}

INITIAL {
	T  = celsius + 273.15
	q10= 1.0^( (celsius-35.0) / 10.0(K) )
	rates()
	xs= xinf
	ys= yinf
}

PROCEDURE rates() { LOCAL a, b  
	a = q10*exp( (1.0e-3)*  zettax*(v-vhalfx)*FARADAY/(R*T) )
	b = q10*exp( (1.0e-3)* -zettax*(v-vhalfx)*FARADAY/(R*T) )
	xinf = a / ( a + b )
	xtau = taox

	a = q10*Ky*exp( (1.0e-3)*  zettay*     gammay *(v-vhalfy)*FARADAY/(R*T) )
	b = q10*Ky*exp( (1.0e-3)* -zettay*(1.0-gammay)*(v-vhalfy)*FARADAY/(R*T) )
	yinf = a   / ( a + b )
	ytau = 1.0 / ( a + b ) + taoy
}
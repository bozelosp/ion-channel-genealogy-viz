NEURON {
	SUFFIX KdBG
	USEION k WRITE ik
	RANGE  gbar,ik
	GLOBAL xtau, ytau, xinf, yinf
}

UNITS {
	(S)	= (siemens)
	(mA)	= (milliamp)
	(mV)	= (millivolt)
	FARADAY	= 96480 (coulombs)
	R	= 8.314  (joule/degC)
}

PARAMETER {
	v		(mV)
	gbar	= 1.0e-3	(S/cm2)
	celsius	= 25	(degC)
	Kx = 1   (1/ms)
	Ky	=   1e-3	(1/ms)
	zettax	=  2.5		(1)
	zettay	=  -1.5		(1)
	vhalfx	= -48.0		(mV)
	vhalfy	= -90.0		(mV)
	taux	=   1		(ms)
	tauy	=   100		(ms)
	q10	= 1.0	(1)    
	FRT = 39 (coulombs/joule) 
}

ASSIGNED {
	ik     	(mA/cm2)
	xtau    (ms)
	ytau    (ms)
	xinf	(1)
	yinf	(1)
}

STATE { xs ys }

BREAKPOINT { 
	SOLVE states METHOD cnexp
	ik= gbar * xs * ys * ( v + 90.0 ) 
}

DERIVATIVE states {
	rates(v)
	xs'= (xinf- xs)/ xtau	
	ys'= (yinf- ys)/ ytau
}

INITIAL {
	rates(v)
	xs= xinf
	ys= yinf
}

PROCEDURE rates(v (mV)) { LOCAL a, b, T, qt
	T = celsius + 273.15  
	qt = q10 ^( (celsius-35.0) / 10.0(K) )
	a = qt*Kx*exp( (1.0e-3)*  zettax*(v-vhalfx)*FRT )
	b = qt*Kx*exp( (1.0e-3)* -zettax*(v-vhalfx)*FRT )
	xinf = a / ( a + b )
	xtau = 1 /(a + b)+ taux

	a = qt*Ky*exp( (1.0e-3)*  zettay* (v-vhalfy)*FRT )
	b = qt*Ky*exp( (1.0e-3)* -zettay* (v-vhalfy)*FRT )
	yinf = a   / ( a + b )
	ytau = 1.0 / ( a + b ) + tauy
}
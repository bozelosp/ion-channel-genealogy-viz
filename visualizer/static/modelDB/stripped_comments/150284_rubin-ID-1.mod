INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}

NEURON {
	SUFFIX rubin
	USEION ca READ cai 
	USEION cal READ cali 
	RANGE P, V, A, B, D, W
}

UNITS {
	(molar) = (1/liter)	
	(mM)	= (millimolar)
	(uM)	= (micromolar)
}

PARAMETER {
	scale = 1	
	pHC = 100		(uM)	
	pHN = 100
	aHC = 28 		(uM)	
	aHN = 3
	
	vtheta = 35		(uM)	
	dtheta = 0.36 			
	btheta = 0.375 			
	
	sigmav = -0.05	(uM)
	sigmad = -0.01
	sigmab = -0.02
	
	tp = 50 	(ms)
	ta = 2		(ms)
	tv = 10		(ms)
	td = 250	(ms)
	tb = 20		(ms)
	tw = 500 	(ms)

	kp = -0.1
	kd = -0.002

	av = 2
	ad = 1.0
	ab = 5.0
	aw = 0.8
	
	bw = 0.6
	
	cp = 1
	cd = 500
	
	p = 0.3
	d = 0.01	
}

STATE {
	A
	V
	B
	D
	P
	W
	
	cai		(mM) 
	cali		(mM) 
}

INITIAL {
	W = 0
	scale_it()
}

ASSIGNED {
	asig
	vsig
	psig
	dsig
	bsig

	w1
	w2
	catot (mM)
}
	
BREAKPOINT {
	catot = cai + cali
	settables(catot)
	SOLVE state METHOD derivimplicit
}

DERIVATIVE state { 

	A' = (asig - A)/ta
	V' = (vsig - V)/tv
	bsig = ab/( 1 + exp((A-btheta)/sigmab) ) 
	B' = (bsig - B - cd*B*V)/tb
	dsig = ad/( 1 + exp((B-dtheta)/sigmad) ) 
	D' = (dsig - D)/td
	P' = (psig - cp*A*P)/tp

	w1 = aw / (1 + exp((P-p)/kp) )
	w2 = bw / (1 + exp((D-d)/kd) )
	W' = (w1 - w2 - W) / tw
}

PROCEDURE settables( cai (mM) ) {
	TABLE asig, vsig, psig DEPEND aHC, aHN, pHC, pHN, vtheta, sigmav, av
		FROM 0.0005 TO 0.1 WITH 200

		asig = ( (cai*(1000)/aHC)^aHN )/ (1 + (cai*(1000)/aHC)^aHN ) 
		vsig = av/( 1 + exp((cai*(1000)-vtheta)/sigmav) ) 
		psig = 10 * ( (cai*(1000)/pHC)^pHN ) / (1 + (cai*(1000)/pHC)^pHN ) 
}

PROCEDURE scale_it() {
	pHC = 4 (uM) * scale 
	aHC = 0.6 (uM) * scale
	
	vtheta = 2 (uM) * scale
	dtheta = 2.6 * scale
	btheta = 0.55 * scale
}
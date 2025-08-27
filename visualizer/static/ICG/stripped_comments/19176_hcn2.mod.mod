NEURON {
	SUFFIX hcn2
	NONSPECIFIC_CURRENT i
	RANGE i, ehcn, g, gbar
	USEION a READ ai  VALENCE 0
	GLOBAL eh
}

UNITS {
	(mV) = (millivolt)
	(molar) = (1/liter)
	(mM) = (millimolar)
}

PARAMETER {
	gbar = 1	(millimho/cm2)
	
	a0 = 0.0015	(/ms)		
	b0 = 0.02	(/ms)
	ah = -135.7    (mV)
	bh = -99.7    (mV)		
	ac = -0.155     (/mV)
	bc = 0.144     (/mV)
	aa0 = 0.0067	(/ms)		
	ba0 = 0.014	(/ms)
	aah = -142.28   (mV)
	bah = -83.5   (mV)
	aac = -0.075    (/mV)
	bac = 0.144 (/mV)
	kon = 3085.7    (/mM-ms)		
	koff = 0.000044857  (/ms)
	b  = 80
	bf = 8.94
	ai	(millimolar)      
	gca = 1			
	shift = 0	(mV)		
	q10v = 4                        
	q10a = 1.5			
	celsius         (degC)
}

ASSIGNED {
	eh (mV)
	v	(mV)
	g	(millimho/cm2)
	i	(milliamp/cm2)
	alpha	(/ms)
	beta    (/ms)
	alphaa	(/ms)
	betaa	(/ms)
}

STATE { c cac o cao }

INITIAL { SOLVE kin STEADYSTATE sparse }

BREAKPOINT {
	SOLVE kin METHOD sparse
	g = gbar*(o + cao*gca)
	i = g*(v-eh)*(1e-3)
}

KINETIC kin {
	LOCAL qa
	qa = q10a^((celsius-22 (degC))/10 (degC))
	rates(v)
	~ c <-> o       (alpha, beta)
	~ c <-> cac  (kon*qa*ai/bf,koff*qa*b/bf)
	~ o <-> cao      (kon*qa*ai, koff*qa)
	~ cac <-> cao      (alphaa, betaa)
	CONSERVE c + cac + o + cao = 1
}

PROCEDURE rates(v(mV)) {
	LOCAL qv
	qv = q10v^((celsius-22 (degC))/10 (degC))
	alpha = a0*qv / (1 + exp(-(v-ah-shift)*ac))
	beta = b0*qv / (1 + exp(-(v-bh-shift)*bc))
	alphaa = aa0*qv / (1 + exp(-(v-aah-shift)*aac))
	betaa = ba0*qv / (1 + exp(-(v-bah-shift)*bac))
}
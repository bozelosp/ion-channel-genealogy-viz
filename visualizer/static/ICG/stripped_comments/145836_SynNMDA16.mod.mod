INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}

NEURON {
	POINT_PROCESS NMDA16
	USEION na READ nao
	RANGE T_max, T, tau, tRel, Erev, synon
	RANGE R,RA,RA2,RA2d1,RA2d2,RA2f,RA2s,O,OMg,RMg,RAMg,RA2Mg,RA2d1Mg,RA2d2Mg,RA2fMg,RA2sMg
	RANGE g, kd1F, kd1B, kd2F, kd2B, csi
	GLOBAL Kcs, kP, kNo, kNi, kMgF, kMgB, ksF,	ksB, kfF, kfB
	NONSPECIFIC_CURRENT i
	THREADSAFE
}

UNITS {
	(nA) = (nanoamp)
	(mV) = (millivolt)
	(pS) = (picosiemens)
	(uS) = (microsiemens)
	(umho) = (micromho)
	(mM) = (milli/liter)
	(uM) = (micro/liter)
}

PARAMETER {

	Erev	= 0    	(mV)	
	gmax	= 50  	(pS)	
	Mg		= 1  	(mM)	
	
	
	tau  = .3 (ms) <1e-9,1e9>
	T_max = 1.5 (mM)		


	kon  = 2.83		(/ms /mM)
	koff = 38.1e-3	(/ms)
	
	ksF0 = 48e-3	(/ms)
	ksB0 = 230e-3	(/ms)
	kfF0 = 2836e-3	(/ms)
	kfB0 = 175e-3	(/ms)
	Vdep = 175		(mV)	
	V0   = -100		(mV)	
	
	
	kd1F = 1e-3		(/ms)
	kd1B = 1e-3		(/ms)
	kd2F = 1e-3		(/ms)
	kd2B = 1e-3		(/ms)

	Kna  = 34.4		(mM)
	Kcs0 = 0.27		(mM)
	a    = -21		(mV)
	kP0  = 1.10e3	(/ms /mM)
	b    = -55		(mV)
	kNo0 = 1.10e2	(/ms)
	c    = 52.7		(mV)
	kNi0 = 61.8e-3	(/ms)
	d	 = -50		(mV)
	csi  = 148		(mM)	
	
	
}

ASSIGNED {
	v		(mV)	
	i 		(nA)	
	g 		(uS)	
	
	T		(mM)	
	tRel	(ms)	
	synon			
	w				
	
	
	ksF		(/ms)
	ksB		(/ms)
	kfF		(/ms)
    kfB		(/ms)
	
	kMgF	(/ms /mM)
	kMgB	(/ms)
	
	Kcs		(mM)
	kP		(/ms /mM)
	kNo		(/ms)
	kNi		(/ms)
	nao		(mM)	
	
	
}

STATE {
	
	R
	RA
	RA2
	RA2d1
	RA2d2
	RA2f
	RA2s
	O
	OMg
	RMg
	RAMg
	RA2Mg
	RA2d1Mg
	RA2d2Mg
	RA2fMg
	RA2sMg
}

INITIAL {
	rates(v)
	T = 0
	synon = 0
	tRel = 0
	
	R = 1
}

BREAKPOINT {
	SOLVE kstates METHOD sparse

	g = w * gmax * O
	i = g * (v - Erev)
}

KINETIC kstates {
	rates(v,t)
	
    
    
	
    
    
	~ R 	 <-> RA			((2*kon*T),koff)
	~ RA 	 <-> RA2		((kon*T),(2*koff))
	~ RA2 	 <-> RA2d1		(kd1F,kd1B)
	~ RA2 	 <-> RA2d2		(kd2F,kd2B)
	~ RA2 	 <-> RA2f		(kfF,kfB)
	~ RA2 	 <-> RA2s 		(ksF,ksB)
	~ RA2f 	 <-> O 			(ksF,ksB)
	~ RA2s 	 <-> O			(kfF,kfB)
	~ O 	 <-> OMg 		((kMgF*Mg),kMgB)
	~ OMg 	 <-> RA2fMg		(ksB,ksF)
	~ OMg 	 <-> RA2sMg		(kfB,kfF)
	~ RA2fMg <-> RA2Mg		(kfB,kfF)
	~ RA2sMg <-> RA2Mg		(ksB,ksF)
	~ RA2Mg  <-> RA2d1Mg	(kd1B,kd1F)
	~ RA2Mg  <-> RA2d2Mg	(kd2B,kd2F)
	~ RA2Mg  <-> RAMg		((2*koff),(kon*T))
	~ RAMg	 <-> RMg		(koff,(2*kon*T))

	CONSERVE R+RA+RA2+RA2d1+RA2d2+RA2f+RA2s+O+OMg+RMg+RAMg+RA2Mg+RA2d1Mg+RA2d2Mg+RA2fMg+RA2sMg = 1
}

NET_RECEIVE(weight) {
	if (flag == 0) {
		tRel = t	
		synon = 1	
					
		w = weight
	}
}

PROCEDURE rates(v (mV), t(ms)) {
	T = T_max * (t - tRel) / tau * exp(1 - (t - tRel) / tau) * synon

	Kcs  = Kcs0 * exp(v/a)
	kP	 = kP0  * exp(v/b)
	kNo  = kNo0 * exp(v/c)
	kNi  = kNi0 * exp(v/d)
	 
	kMgF = kP / ((1 + nao/Kna) * (1 + nao/Kna + csi/Kcs))
	kMgB = kNo / (1 + nao/Kna)^2 + kNi
	

	ksF = ksF0 * exp((v - V0) / Vdep)
	ksB = ksB0 * exp((v - V0) / Vdep * (-1))
	kfF = kfF0 * exp((v - V0) / Vdep)
	kfB = kfB0 * exp((v - V0) / Vdep * (-1))
}
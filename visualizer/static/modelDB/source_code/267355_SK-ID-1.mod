TITLE Calcium activated Potassium channel (SK)  

: Author: Chitaranjan Mahapatra (chitaranjan@iitb.ac.in)
: Computational Neurophysiology Lab
: Indian Institute of Technology Bombay, India 

: For details refer: 
: Mahapatra C, Brain KL, Manchanda R, A biophysically constrained computational model of the action potential 
: of mouse urinary bladder smooth muscle. PLOS One (2018) 
:Kinetics were fit to data from Filipis et al. 2022

:SK channel ( function of Calcium only)

NEURON {
	SUFFIX sk
	USEION k READ ek WRITE ik
        USEION tca READ tcai
        USEION nca READ ncai
        RANGE  gbar,gkahp,ik, ninf,taun,g
        GLOBAL Cq10,mt, a0, b0, cahco, carco
}

UNITS {
	(mA) = (milliamp)
	(mV) = (millivolt)
	(molar) = (1/liter)
	(mM) = (millimolar)
	(pS) = (picosiemens)
	(um) = (micron)
}

PARAMETER {
	gbar = 0.01	(S/cm2)
       : n = 4
        cai = 50.e-6	(mM)
        a0 = 1.3e4	(1/ms-mM)	
        b0 = 0.06	(1/ms)				
	    celsius = 37(degC)
	Cq10 = 2
	mt = 0.2
	
	tcai (mM)
	ncai (mM)
	cahco=0.0
	carco=1
	
}

STATE {	n }

ASSIGNED {
	ik	(mA/cm2)
        g	(S/cm2)
        ninf
        taun	(ms)
	a	(1/ms)
        v	(mV)
        ek	(mV)
}

BREAKPOINT {
	SOLVE state METHOD cnexp
	g = gbar*n
	ik = g*(v-ek)
}

INITIAL {
	rate(tcai, ncai)
	n=ninf
}

DERIVATIVE state {
	rate(tcai,ncai)
	n' = (ninf - n)/taun
}


PROCEDURE rate(tcai (mM), ncai(mM)) {
	LOCAL q10
	q10 = Cq10^((celsius - 22 (degC))/10 (degC) )
	a = a0*(tcai*carco+ncai*cahco)/10
		if (a < 0.05) {
	       a = 0 
	} 
	taun = b0
	ninf = a
}



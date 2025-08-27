NEURON {
	SUFFIX cldif_CA3
	USEION cl READ icl WRITE cli VALENCE -1
	USEION hco3 READ hco3i, hco3o VALENCE -1
	GLOBAL vrat		
	RANGE tau, cli0, clo0, egaba, delta_egaba, init_egaba, ehco3_help, ecl_help
}

DEFINE Nannuli 4

UNITS {
	(molar) = (1/liter)
	(mM) = (millimolar)
	(um) = (micron)
	(mA) = (milliamp)
	(mV)    = (millivolt)
	FARADAY = (faraday) (10000 coulomb)
	PI = (pi) (1)
	F = (faraday) (coulombs)
	R = (k-mole)  (joule/degC)
}

PARAMETER {
	DCl = 2 (um2/ms) 
	tau = 174000 (ms) 
	cli0 = 8 (mM)
	clo0 = 133.5 (mM)
	hco3i0 = 16	(mM)
	hco3o0 = 26	(mM)
	P_help = 0.18
	celsius = 37    (degC)

}

ASSIGNED {
	diam 	(um)
	icl 	(mA/cm2)
	cli 	(mM)
	hco3i	(mM)
	hco3o	(mM)
	vrat[Nannuli]	
			
			
	egaba 	(mV)
	ehco3_help 	(mV)
	ecl_help	(mV)
	init_egaba  (mV)
	delta_egaba (mV)
}

STATE {
	
	
	cl[Nannuli]	(mM) <1e-10>
}


BREAKPOINT { 
		SOLVE state METHOD sparse
		ecl_help = log(cli/clo0)*(1000)*(celsius + 273.15)*R/F
		egaba = P_help*ehco3_help + (1-P_help)*ecl_help
		delta_egaba = egaba - init_egaba
}

LOCAL factors_done

INITIAL {
	if (factors_done == 0) {  	
		factors_done = 1	
		factors()		
	}
	cli = cli0
	hco3i = hco3i0
	hco3o = hco3o0
	FROM i=0 TO Nannuli-1 {
		cl[i] = cli
	}
	ehco3_help = log(hco3i/hco3o)*(1000)*(celsius + 273.15)*R/F
	ecl_help = log(cli/clo0)*(1000)*(celsius + 273.15)*R/F
	egaba = P_help*ehco3_help + (1-P_help)*ecl_help
	init_egaba = egaba
	delta_egaba = egaba - init_egaba 
}

LOCAL frat[Nannuli]	

PROCEDURE factors() {
	LOCAL r, dr2
	r = 1/2			
	dr2 = r/(Nannuli-1)/2	
				
	vrat[0] = 0
	frat[0] = 2*r		
	FROM i=0 TO Nannuli-2 {
		vrat[i] = vrat[i] + PI*(r-dr2/2)*2*dr2	
		r = r - dr2
		frat[i+1] = 2*PI*r/(2*dr2)	
						
		r = r - dr2
		vrat[i+1] = PI*(r+dr2/2)*2*dr2	
	}
}

KINETIC state {
	COMPARTMENT i, diam*diam*vrat[i] {cl}
	LONGITUDINAL_DIFFUSION i, DCl*diam*diam*vrat[i] {cl}
	~ cl[0] << ((icl*PI*diam/FARADAY) + (diam*diam*vrat[0]*(cli0 - cl[0])/tau)) 
	FROM i=0 TO Nannuli-2 {
		~ cl[i] <-> cl[i+1]	(DCl*frat[i+1], DCl*frat[i+1])
	}
	cli = cl[0]
}
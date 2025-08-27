NEURON {
	SUFFIX cldif_CA3_NKCC1_HCO3
	USEION cl READ icl WRITE cli VALENCE -1 
	USEION hco3 READ ihco3 WRITE hco3i VALENCE -1
	GLOBAL vrat		
	RANGE tau, cli0, clo0, hco3i0, hco3o0, egaba, delta_egaba, init_egaba, ehco3_help, ecl_help 
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
	tau_NKCC1 = 174000 (ms)   
	tau_passive = 321000 (ms) 
        tau_hco3 = 1000 (ms) 
	cli0 = 50 (mM) 
	cli_Start = 10 (mM) 
	clo0 = 133.5 (mM) 
	hco3i0 = 16	(mM) 
	hco3o0 = 26	(mM) 
	hco3i_Start = 16 (mM) 
	celsius = 31    (degC)

}

ASSIGNED {
	diam 	(um)
	icl 	(mA/cm2) 
        ihco3 	(mA/cm2) 
	cli 	(mM) 
	hco3i	(mM) 
	hco3o	(mM) 
	vrat[Nannuli]	
			
			
	ehco3_help 	(mV)
	ecl_help	(mV)
	ActPump   
}

STATE {
	
	
	cl[Nannuli]	(mM) <1e-10>
        hco3[Nannuli]	(mM) <1e-10>
}


BREAKPOINT {
		SOLVE state METHOD sparse
		ecl_help = log(cli/clo0)*(1000)*(celsius + 273.15)*R/F
                ehco3_help = log(hco3i/hco3o0)*(1000)*(celsius + 273.15)*R/F
}

LOCAL factors_done

INITIAL {
	if (factors_done == 0) {  	
		factors_done = 1	
		factors()		
	}
	cli = cli_Start
	hco3i = hco3i0
	hco3o = hco3o0
	FROM i=0 TO Nannuli-1 { 
		cl[i] = cli
	}
        FROM i=0 TO Nannuli-1 { 
		hco3[i] = hco3i
	}
	ehco3_help = log(hco3i/hco3o)*(1000)*(celsius + 273.15)*R/F 
	ecl_help = log(cli/clo0)*(1000)*(celsius + 273.15)*R/F
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
    if (cli0 >= cl[0]) { 
		  ActPump = 1
		}
		else {     
		  ActPump = 0
		}

  	COMPARTMENT i, diam*diam*vrat[i] {cl}
		LONGITUDINAL_DIFFUSION i, DCl*diam*diam*vrat[i] {cl}
				~ cl[0] << ((icl*PI*diam/FARADAY) + ActPump*(diam*diam*vrat[0]*(cli0 - cl[0])/tau_NKCC1) + (diam*diam*vrat[0]*(cli0 - cl[0])/tau_passive)) 
	 	FROM i=0 TO Nannuli-2 {
		~ cl[i] <-> cl[i+1]	(DCl*frat[i+1], DCl*frat[i+1])
                }
	        cli = cl[0]

        COMPARTMENT i, diam*diam*vrat[i] {hco3}
		LONGITUDINAL_DIFFUSION i, DCl*diam*diam*vrat[i] {hco3}
				~ hco3[0] << ((ihco3*PI*diam/FARADAY)  + (diam*diam*vrat[0]*(hco3i0 - hco3[0])/tau_hco3)) 
	 	FROM i=0 TO Nannuli-2 {
		~ hco3[i] <-> hco3[i+1]	(DCl*frat[i+1], DCl*frat[i+1])
                }
	        hco3i = hco3[0]
}
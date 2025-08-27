NEURON {
		SUFFIX NaCaPump									
		USEION ca READ cao, cai WRITE ica				
		USEION na READ nao, nai WRITE ina				
		RANGE  inca, DFout, DFin, S, KNaCa, DNaCa		
	}


	UNITS {
		(mA) = (milliamp)
		(mV) = (millivolt)
		(molar) = (1/liter)
		(mM) = (millimolar)
		FARADAY = 96500 (coulombs)
		R = 8.314 (joule/degC)
	}


	PARAMETER {
		KNaCa22 = 1.27324E-06 (mA/cm2/mM4)    <0,1e6> 
		Q10NaCa = 2.20				
		Q10TempA = 22.85	(degC)		
		Q10TempB = 10	(degC)
		r=3
		gamma=0.5
		DNaCa=0.0036 (/mM4)
	}


	ASSIGNED {
		
		
		celsius (degC)
		v (mV)
		cai (mM)
		cao (mM)
		ica (mA/cm2)
		ina (mA/cm2)
		nao (mM)
		nai (mM)
		
		
		inca (mA/cm2)
		S
		DFin (mM4)
		DFout (mM4)
		temp (degC)
		KNaCa (mA/cm2/mM4) 
	}


	BREAKPOINT {

		temp = celsius +273.15

		S=1.0+DNaCa*(cai*nao*nao*nao+cao*nai*nai*nai)
		
		DFin=nai*nai*nai*cao*exp(((r-2)*gamma*v*FARADAY)/((1000)*R*temp))
		
		DFout=nao*nao*nao*cai*exp(((r-2)*(gamma-1)*v*FARADAY)/((1000)*R*temp))

		inca=KNaCa*((DFin-DFout)/S)
		
		
		ina = 3*inca
		ica = -2*inca
	}
	

	INITIAL {
		KNaCa = KNaCa22*Q10NaCa^((Q10TempA-celsius)/Q10TempB)
	}
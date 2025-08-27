NEURON {
		SUFFIX leakSchild							
		USEION na READ ena WRITE ina		
		USEION ca READ cai, cao WRITE ica	
		RANGE i, ina, ica, gbna, gbca		
	}


	UNITS {
		(S)=(siemens)
		(mV)=(millivolt)
		(mA)=(milliamp)
		FARADAY = 96500 (coulombs)
		(molar) = (1/liter)
		(mM)    = (millimolar)
	}


	PARAMETER {
		gbna=1.85681E-05 (S/cm2) <0, 1e9>
		gbca=3.00626E-06 (S/cm2) <0, 1e9>
		
		R=8.314 (joule/degC)
		z=2 
		ecaoffset=78.7 (mV)
	}


	ASSIGNED {

		
		v (mV)
		ena (mV)
		ica (mA/cm2)
		ina (mA/cm2)
		celsius (degC)
		cai(mM)
		cao (mM)
		
		
		ecaleak	(mV)
	}


	BREAKPOINT { 
		
		
		ecaleak=(1000)*(R*(celsius+273.15)/z/FARADAY*log(cao/cai))-ecaoffset 
		
		ina = gbna*(v - ena) 
		ica = gbca*(v - ecaleak) 
	}
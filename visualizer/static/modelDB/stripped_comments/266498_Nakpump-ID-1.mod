NEURON {
	SUFFIX nakpump
	USEION na READ nai WRITE ina
	USEION k READ ko WRITE ik
	RANGE ik, ina
	RANGE nain, naout, kin, kout, smalla, b1, b2, celsiusT, kvotqt
}

UNITS {
        (molar) = (1/liter)                     
        (mM)    = (millimolar)             	
	
	(uA) = (microamp)
	(mA) = (milliamp)
	(mV) = (millivolt)
	(um) = (micron)

}

PARAMETER {


	nain = 11.4 	(mM)
	kout = 5.6 	(mM)

	smalla = 0  	(mA/cm2)
	b1 = 1  	(mM)

	
	
	kvotqt

	celsiusT = 32

}

ASSIGNED {

	ina (mA/cm2)
	ik (mA/cm2)

	nai (mM)
	ko (mM)
}

STATE {
	inapump  (mA/cm2)
	ikpump (mA/cm2)	

}


INITIAL {

	inapump = 0
	ikpump = 0

}

BREAKPOINT {


	kvotqt = 1^((celsius-22)/10)

	ikpump = smalla/((1+b1/ko)^2) * (1.62/(1+(6.7(mM)/(nai+8(mM)))^3) + 1.0/(1+(67.6(mM)/(nai+8(mM)))^3))
	ikpump = ikpump*kvotqt
	inapump = -1.5*ikpump

	ina = inapump		
	ik = ikpump		

}
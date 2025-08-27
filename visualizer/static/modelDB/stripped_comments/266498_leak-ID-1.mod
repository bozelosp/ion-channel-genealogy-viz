NEURON {
	SUFFIX leak
	USEION k READ ek WRITE ik
	USEION na READ ena WRITE ina
	USEION ca READ eca WRITE ica
	RANGE gkleak, gnaleak, gk, qna, ik, ina, ica
}

UNITS { 
	(mV) = (millivolt)  (mA) = (milliamp)
	(um) = (micron)
	PI		= (pi) (1)
	FARADAY		= 96485.309 (coul)
}

PARAMETER {


	gkleak	= 0 (mho/cm2)
	gnaleak	= 0 (mho/cm2)
	gcaleak = 0 (mho/cm2)
	
	
}

ASSIGNED {
	v (mV)
	ik (mA/cm2)
	ek (mV)
	ina (mA/cm2)
	ena (mV)
	ica (mA/cm2)
	eca (mV)

}

BREAKPOINT {
	ik = gkleak*(v-ek)
	ina = gnaleak*(v-ena)
	ica = gcaleak*(v-eca)
}


INITIAL {
	ik = gkleak*(v-ek)
	ina = gnaleak*(v-ena)
	ica = gcaleak*(v-eca)
}
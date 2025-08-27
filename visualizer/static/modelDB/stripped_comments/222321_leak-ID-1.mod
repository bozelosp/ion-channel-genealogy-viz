NEURON {
	SUFFIX leak
	USEION k READ ek WRITE ik
	USEION na READ ena WRITE ina
        NONSPECIFIC_CURRENT ifix
	RANGE gk, ik, gna, ina, gfix, ifix
        RANGE qk, qna
}

UNITS { 
	(mV) = (millivolt)
        (mA) = (milliamp)
        PI		= (pi) (1)
	FARADAY		= 96485.309 (coul/mole)
}

PARAMETER {
	gk	= 1e-5 (mho/cm2)  
	gna	= 1e-5 (mho/cm2)  
	gfix	= 5e-4 (mho/cm2)  
}

ASSIGNED {
	v (mV)
        v_init (mV)
	ik (mA/cm2)
	ek (mV)
	ina (mA/cm2)
	ena (mV)
	ifix (mA/cm2) 
        diam (um)
}

BREAKPOINT {
	ik = gk*(v-ek)
	ina = gna*(v-ena)
	
	ifix = gfix*(v+70)
        SOLVE integrate METHOD sparse
}

STATE { qk qna }

INITIAL {
	ik = gk*(v-ek)
	ina = gna*(v-ena)
	qk = 0
	qna = 0
}

KINETIC integrate {
	COMPARTMENT diam*diam*PI/4 { qna qk }
	~ qna << ((-ina*diam)*PI*(1e4)/FARADAY )
	~ qk  << ((-ik*diam)*PI*(1e4)/FARADAY )
}
UNITS {
      (mv) = (millivolt)
      (mA) = (milliamp)
}

NEURON {
       SUFFIX NaL
       USEION na READ ena,nai WRITE ina
       RANGE gna,inaL
       GLOBAL activate_Q10,gmaxQ10,gmax_k,temp1,temp2,tempb
}

PARAMETER {
        v (mV)
	gna = 0.81e-5 (mho/cm2)
	inaL = 0.0 (mA/cm2)
	ena
	nai
	celsius

	activate_Q10 = 1
	gmaxQ10 = 1.5
	temp1 = 25.0 (degC)
	temp2 = 35.0 (degC)
	tempb = 23.0 (degC)	
}

ASSIGNED { 
        ina (mA/cm2)
	gmax_k
}

BREAKPOINT {
	   ina	= gna*gmax_k*(v-ena)
	   inaL = ina
}
UNITSOFF

INITIAL {
	LOCAL ktemp,ktempb,ktemp1,ktemp2
	if (activate_Q10>0) {
	  ktemp  = celsius+273.0
	  ktempb = tempb+273.0
	  ktemp1 = temp1+273.0
	  ktemp2 = temp2+273.0
	  gmax_k = exp( log(gmaxQ10)*((1/ktempb)-(1/ktemp))/((1/ktemp1)-(1/ktemp2)) )
	}else{
	  gmax_k = 1.0
	}
}	
UNITSON
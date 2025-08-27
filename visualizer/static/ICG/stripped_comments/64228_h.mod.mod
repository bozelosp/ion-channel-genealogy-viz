NEURON 
{
	SUFFIX h
	
	
	NONSPECIFIC_CURRENT i
	RANGE  gh, ghbar,  Vhalfh, Sh
	RANGE  aoh
	GLOBAL eh
	

}

UNITS
{
	(mA) = (milliamp)
	(mV) = (millivolt)
	(mS) = (millimho)
}

PARAMETER
{
	
	ghbar = 1.41 (mS/cm2) <0,1e9>
	
        
       
        Vhalfh=-82  (mV)
       
        Sh = -5.33    (mV)
       
	aoh  = 1.64 (/s)
	
	
	

}

STATE
{
	
	nh
	
}

ASSIGNED
{
	eh (mV)
	v (mV)
	
	i (mA/cm2)
	
	infh
	tauh   (ms)
	
	gh (mho/cm2)
	


}

INITIAL
{
	rate(v)
	
	nh  = infh
	
}

BREAKPOINT
{
	SOLVE states METHOD cnexp
	
	gh  = (0.001)*ghbar*nh
	i  = gh*(v - eh)
	
	
}

DERIVATIVE states
{
	rate(v)
	
	nh'  = (infh  - nh )/tauh
	
}




FUNCTION alphah(v(mV)) (/ms)
{ 
	alphah = 0.001*aoh*exp( (v - Vhalfh)/(2*Sh)   )

}

FUNCTION betah(v(mV)) (/ms)
{ 
	betah =0.001*aoh*exp( - ( v-Vhalfh)/(2*Sh) )

}




PROCEDURE rate(v (mV))
{
        LOCAL  ah, bh

	
	
	ah = alphah(v)
	bh = betah(v)
	tauh = 1/(ah + bh)
	infh = ah*tauh
	
	

}
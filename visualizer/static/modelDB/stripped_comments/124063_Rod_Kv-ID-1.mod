NEURON 
{
	SUFFIX Kv
	
	USEION Kv WRITE iKv VALENCE 1
		
        RANGE  gKvbar, eKv
	
}

UNITS
{
	(mA) = (milliamp)
	(mV) = (millivolt)
	(mS) = (millimho)
	
}

PARAMETER
{
        
        gKvbar = 10 (mS/cm2) <0,1e9> 
        eKv = -80 (mV) 

}

STATE
{
	mKv
	
}

ASSIGNED
{
	v (mV)
	
	iKv (mA/cm2)
              
          
	infmKv
	taumKv  (ms)
	  
	gKv (mho/cm2)
	
}

INITIAL
{
	rate(v)
	mKv = infmKv
	
}




BREAKPOINT
{
	SOLVE states METHOD cnexp
	gKv = (0.001)*gKvbar*mKv*mKv*mKv*mKv
	iKv = gKv*(v - eKv)
}

DERIVATIVE states
{
	rate(v)
	mKv' = (infmKv - mKv)/taumKv
	
}


UNITSOFF

FUNCTION alphamKv(v(mV)) (/ms)
{ 
	alphamKv = 0.005*(20-v)/( exp( (20-v)/22) -1 )
 
}

FUNCTION  betamKv (v(mV)) (/ms)
{
	
	betamKv = 1/16*exp (- v/80 )
}



UNITSON

PROCEDURE rate(v (mV))
{
        LOCAL am, bm

	
	am = alphamKv(v)
	bm = betamKv(v)
	taumKv = 1/(am + bm)
	infmKv = am*taumKv

}
NEURON 
{
	SUFFIX Kv
	
	
	USEION k READ ek WRITE ik
	
        RANGE  gKvbar
	
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
        

}

STATE
{
	mKv
	
}

ASSIGNED
{
	v (mV)
	ek (mV)
	ik (mA/cm2)
              
          
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
	ik = gKv*(v - ek)
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
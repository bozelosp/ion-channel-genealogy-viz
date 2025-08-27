NEURON 
{
	SUFFIX hyper
	
	NONSPECIFIC_CURRENT i
	GLOBAL eh
	RANGE ghyperbar, ghyper, ehyper 
         

}

UNITS
{
	(mA) = (milliamp)
	(mV) = (millivolt)
	(mS) = (millimho)
}

PARAMETER
{
       
   
       
        ghyperbar = 3.5 (mS/cm2) <0,1e9>
        ehyper  = -32.5 (mV)
       
       
        

}

STATE
{
	nhyper
	
}

ASSIGNED
{
	v (mV)
	eh (mV)
	i (mA/cm2)
	
	infhyper
	tauhyper   (ms)
	
	ghyper (mho/cm2)

}

INITIAL
{
	rate(v)
	nhyper  = infhyper
}




BREAKPOINT
{
	SOLVE states METHOD cnexp
	ghyper  = (0.001)*ghyperbar*(1-(1+3*nhyper)*(1-nhyper)^3)
	i  = ghyper*(v - eh)
	
	
	
	
}

DERIVATIVE states
{
	rate(v)
	nhyper'  = (infhyper  - nhyper )/tauhyper

}



FUNCTION alphah(v(mV)) (/ms)
{ 
	alphah = 0.001*18/( exp  (  ( v+88)/12 ) + 1 )
}


FUNCTION betah(v(mV)) (/ms)
{ 
	betah = 0.001*18/( exp  ( - ( v+18)/19 ) + 1 )
}



PROCEDURE rate(v (mV))
{
        LOCAL a, b

	
	a = alphah(v)
	b = betah(v)
	tauhyper = 1/(a + b)
	infhyper = a/(a + b)
	
}
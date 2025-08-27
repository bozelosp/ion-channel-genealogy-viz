NEURON 
{
	SUFFIX CPR
	
	USEION Ca WRITE iCa VALENCE 2
	USEION Cl WRITE iCl  VALENCE 1
	USEION Kca WRITE iKca VALENCE 1
	
	NONSPECIFIC_CURRENT  iCGMP
	
	RANGE gCabar, gCa, eCa, SCa, VhalfCa, aoCa
	
             RANGE gClbar,gCl, eCl, SCl
             RANGE gKcabar,gKca, eKca
	
	
	RANGE gCGMP, eCGMP
	
	
	RANGE FactorCaI
	RANGE mCl, Cas
	
	
	

}

UNITS
{
	(mA) = (milliamp)
	(mV) = (millivolt)
	(mS) = (millimho)
	(mol)= (1)
	(M)  = (mol/liter)
	(uM) = (micro M)
}

PARAMETER
{
       
       
       gCabar = 4.9 (mS/cm2) <0,1e9>
       eCa =  40 (mV)
       aoCa = 0.0031  (/ms)
       VhalfCa=-16.6 (mV)
       SCa =5.7      (mV)   
   
       
        
       
       eCl= -45  (mV)
       gClbar = 6.5 (mS/cm2) <0,1e9>
       SCl = 0.09   (uM)                                                     
       Clh = 0.37 (uM)
       FactorCaI = 0.5    
 
       
       eKca=-80 (mV)
       gKcabar = 0.5 (mS/cm2) 
 
      	
       
        gCGMP= 0   (mS/cm2)
	
        eCGMP=0.8 (mV)
        

}

STATE
{

	nCa
	mKca
	
}

ASSIGNED
{
	v (mV)
	
	iCa (mA/cm2)
             iCl  (mA/cm2)
             iCGMP (mA/cm2) 
             iKca (mA/cm2) 
              
           
           
	infmKca
	taumKca  (ms)
	
	infCa
	tauCa  (ms) 
	
	Cas  (uM)
	mCl
	
	
	mKca1
	
	gKca (mho/cm2)
	
	gCa (mho/cm2)
             gCl (mho/cm2)

}

INITIAL
{
	rate(v)
	nCa = infCa
	mKca= infmKca
}




BREAKPOINT
{
	SOLVE states METHOD cnexp
	gCa = (0.001)*gCabar*nCa
	iCa = gCa*(v - eCa)
	
	UNITSOFF
	
	
	
	
	
	
		Cas =-0.2+FactorCaI * (-iCa) * 1 *  0.5         /(1.6e-19)/  (6.023e23) * 1e-6         *1e14    
	
	
         
                
	
	
	mCl = 1/(1+ exp ( (Clh - Cas)/ SCl  ) ) 
	gCl = (0.001)* gClbar * mCl
	iCl = gCl*(v-eCl)   
	
	mKca1=Cas/(Cas+0.3)
	gKca=(0.001)*gKcabar*mKca*mKca*mKca1
	iKca=gKca*(v-eKca) 
	
	UNITSON
	
	  
	iCGMP = (0.001)*gCGMP*(v-eCGMP)
	
	
	
	
	
}

DERIVATIVE states
{
	rate(v)
	nCa' = (infCa - nCa)/tauCa
	mKca'= (infmKca - mKca ) /taumKca

}


UNITSOFF


FUNCTION alphamKca(v(mV)) (/ms)
{ 
	alphamKca = (0.001)*15*(80-v)/ ( exp( (80-v)/40 ) -1)
	
}

FUNCTION  betamKca (v(mV)) (/ms)
{
	
	betamKca = (0.001)*20*exp (-v/35)
}



UNITSON


FUNCTION alphaCa(v(mV))(/ms)
{ 
	alphaCa = aoCa*exp( (v - VhalfCa)/(2*SCa)   )
}

FUNCTION betaCa(v(mV))(/ms)
{ 
	betaCa = aoCa*exp( - ( v-VhalfCa)/(2*SCa) )
}


PROCEDURE rate(v (mV))
{
        LOCAL a, b

	
	
	a = alphamKca(v)
	b = betamKca(v)
	taumKca = 1/(a + b)
	infmKca = a/(a + b)
	
	
	
	a = alphaCa(v)
	b = betaCa(v)
	tauCa = 1/(a + b)
	infCa = a/(a + b)

}
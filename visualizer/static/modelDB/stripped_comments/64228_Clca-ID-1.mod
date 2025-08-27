NEURON 
{
	SUFFIX Clca
		
	USEION Ca READ Cai VALENCE 2
	
	USEION Cl WRITE iCl  VALENCE 1
	
	RANGE gClbar, eCl, Clh
	
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
       
       eCl= -45  (mV)
       gClbar = 3 (mS/cm2) <0,1e9>
       
       Clh = 0.5 (uM )
       Cai   (mM)
       

}

STATE
{
	mCl
	
}

ASSIGNED
{
	v (mV)
	iCl (mA/cm2)
	
	
        gCl (mho/cm2)

}

INITIAL
{
		
	
}




BREAKPOINT
{       LOCAL Cas
 
	Cas=Cai*1000  	
	mCl = 1/(  1+(Clh/Cas)^4 ) 
	gCl = (0.001)* gClbar * mCl
	iCl = gCl*(v-eCl)   
		
	
}
NEURON 
{
	SUFFIX Kca
	USEION Ca READ Cai VALENCE 2
	USEION Kca WRITE iKca VALENCE 1
	RANGE infmKcaV,taumKcaV,eKca,gKcabar,Cahalf
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
       
 
       
       eKca=-80 (mV)
       gKcabar = 5 (mS/cm2)
       Cai
	Cahalf=0.32 (uM)	 
 

}

STATE
{

	
	mKcaCa
	
}

ASSIGNED
{
	iKca (mA/cm^2)
	v (mV)
           
	
	
	gKca (mho/cm2)

}

INITIAL
{      LOCAL Cas
	
	Cas=Cai*1000 

	
        
        mKcaCa=1/(1+(Cahalf/Cas)^4)
}




BREAKPOINT
{
	LOCAL Cas
	
	Cas=Cai*1000 
	mKcaCa=1/(1+(Cahalf/Cas)^4 )
	
        gKca=(0.001)*gKcabar*mKcaCa^4
	iKca=gKca*(v-eKca) 
	
	
}
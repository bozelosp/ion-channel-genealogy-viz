NEURON 
{
	SUFFIX Kca
	
	
	USEION ca READ cai
        USEION k READ ek WRITE ik
        RANGE infmKcaV,taumKcaV,gKcabar,Cahalf
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
       
 
       
       
       gKcabar = 5 (mS/cm2)
       
	Cahalf=0.32 (uM)	 
 

}

STATE
{

	
	mKcaCa
	
}

ASSIGNED
{
	
	ik (mA/cm^2)
        v (mV)
        ek (mV)
        cai (mM)
           
	
	
	gKca (mho/cm2)

}

INITIAL
{      LOCAL Cas
	
	Cas=cai*1000 

	
        
        mKcaCa=1/(1+(Cahalf/Cas)^4)
}




BREAKPOINT
{
	LOCAL Cas
	
	Cas=cai*1000 
	mKcaCa=1/(1+(Cahalf/Cas)^4 )
	
        gKca=(0.001)*gKcabar*mKcaCa^4
	ik=gKca*(v-ek) 
	
	
}
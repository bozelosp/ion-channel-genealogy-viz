NEURON 
{
	SUFFIX Leak
	

	
	NONSPECIFIC_CURRENT il
	
	
          
	RANGE glbar, el
	
	

	

}

UNITS
{
	(mA) = (milliamp)
	(mV) = (millivolt)
	(mS) = (millimho)
        
}	

PARAMETER
{
       
       
       
        glbar=0.6 (mS/cm2)
       
	
	 el=-55 (mV)

}


ASSIGNED
{
	v (mV)
	gl (mho/cm2)
	il  (mA/cm2)
          
}





BREAKPOINT
{
        gl= 1e-3*glbar 

	il  = gl*(v-el)

	
	
}
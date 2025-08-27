UNITS {
    (mA) = (milliamp)
    (mV) = (millivolt)
	(S)  = (siemens)
}






NEURON 	{
	
	
	
	
	
	
	SUFFIX hcn

	NONSPECIFIC_CURRENT ihcn

	
	
	
	
	
	RANGE ghcn, V12, ihcn, ehcn, Ft

	
	
	
	

}





PARAMETER 	{

        ehcn    = -30 		(mV)
        Ft      = 1			(ms)
    	ghcn	= 0.000085	(S/cm2)
    	V12		= -90		(mV)
}















ASSIGNED {
	
    v 		(mV)
    celsius (degC)
    
    ihcn	(mA/cm2)
    sinf
    stau
}








STATE { shcn }











BREAKPOINT 	{
	
	SOLVE states METHOD cnexp
	
	ihcn = ghcn * shcn * ( v-ehcn )
}





INITIAL 	{
	rates()
	shcn = sinf
}















DERIVATIVE states	{
	rates()
	shcn' = (sinf-shcn)/stau
}







PROCEDURE rates(){
	UNITSOFF 
		sinf = 1.0000/( 1.0000 + exp((v- V12)/8.0000) ) 
		stau = Ft*exp( 0.033*(v+75) )/( 0.013*( 1+exp(0.083*(v+75)) ) )
	UNITSON
}
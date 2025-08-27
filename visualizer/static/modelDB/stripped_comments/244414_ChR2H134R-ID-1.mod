UNITS {

	(mS) = (millisiemens)
	(mV) = (millivolt)
	(mA) = (milliamp)
}


NEURON {

	SUFFIX ChR2
	NONSPECIFIC_CURRENT i

	RANGE i, gmax, Irradiance
    RANGE light_intensity, light_delay, pulse_width

    GLOBAL A, B, C, gamma
	GLOBAL wavelength, hc, wloss, sigma_retinal
	GLOBAL Q10_Gd1, Q10_Gd2, Q10_Gr, Q10_e12dark, Q10_e21dark, Q10_epsilon1, Q10_epsilon2
}



PARAMETER {


	Irradiance      = 0.
	gmax  			= 0.4  		(mS/cm2)	  
	EChR2 			= 0.     	(mV)          

	light_delay     = 100.		(ms)		  
	pulse_width     = 100.		(ms)		  
	light_intensity = 5.					  

    gamma           = 0.1					  
	A 				= 10.6408   (mV)          
	B 				= -14.6408  (mV)		  
	C 				= -42.7671  (mV)          

	wavelength 		= 470		    		  
	hc       		= 1.986446E-25  		  
	wloss    		= 1.3      
	sigma_retinal 	= 12.E-20       		  

	tauChR2  		= 1.3		(ms)          

	temp 			= 22	  (degC)		  
	Q10_Gd1      	= 1.97					  
	Q10_Gd2      	= 1.77					  
	Q10_Gr       	= 2.56					  
	Q10_e12dark  	= 1.1					  
	Q10_e21dark  	= 1.95					  
	Q10_epsilon1 	= 1.46					  
	Q10_epsilon2 	= 2.77					  
}


ASSIGNED {

	v      		(mV)			
	celsius		(degC)          

	i    		(mA/cm2)		

	Gd1    		(1./ms)			
	Gd2    		(1./ms)			
	Gr     		(1./ms)			
	e12    		(1./ms)			
	e21    		(1./ms)			

	epsilon1 					
	epsilon2 					
	F   		(1./ms)			
	S0							
}


STATE {







	O1 O2 C1 p
}


BREAKPOINT {
	SOLVE states METHOD cnexp

	Irradiance = 0.               

    
    
    
    

    if (t < light_delay)                      { Irradiance = 0. }
    else if (t < (light_delay + pulse_width)) { Irradiance = light_intensity }
    else if (t > (light_delay + pulse_width)) { Irradiance = 0. }

	i = (0.001) * (gmax * (A + B * exp(v /C))/v * (O1 + gamma * O2) * (v - EChR2))

	
}


INITIAL {

	rates(v)
	C1 = 1.
	O1 = 0.
	O2 = 0.
	p  = 0.
}


DERIVATIVE states {

    rates(v)
    
	
	
	
	
	

	
	O1' = -(Gd1+e12)           * O1 + e21 * O2 + epsilon1*F*p * C1
	O2' = -(Gd2+e21)           * O2 + e12 * O1 + epsilon2*F*p * (1. - C1 - O1 - O2)
	C1' = -epsilon1*F*p        * C1 + Gd1 * O1 + Gr  * (1. - C1 - O1 - O2)  

	p'  =  (S0 - p) / tauChR2
}


UNITSOFF
PROCEDURE rates(v (mV)) {
    LOCAL e12dark, e21dark, logphi0, Ephoton, flux

	
	e12dark  = 0.011                                 
	e21dark  = 0.008                                 
	epsilon1 = 0.8535
	epsilon2 = 0.14
	Gd1 = 0.075 + 0.043 * tanh( -(v+20.) / 20.)    	 
	Gd2 = 0.05                                  	 
	Gr  = 0.0000434587 * exp(-0.0211539274 * v)    	 

	
	e12dark  = e12dark  * Q10_e12dark^((celsius-temp)/10.)    
	e21dark  = e21dark  * Q10_e21dark^((celsius-temp)/10.)    
	epsilon1 = epsilon1 * Q10_epsilon1^((celsius-temp)/10.)   
	epsilon2 = epsilon2 * Q10_epsilon2^((celsius-temp)/10.)   
	Gd1 	 = Gd1           * Q10_Gd1^((celsius-temp)/10.)	  
	Gd2 	 = Gd2           * Q10_Gd2^((celsius-temp)/10.)	  
	Gr  	 = Gr             * Q10_Gr^((celsius-temp)/10.)   

	if (Irradiance>0) {
		logphi0  = log(1. + Irradiance / 0.024)      
	}
	else {
		logphi0 = 0.	 							 
	}
	e12      = e12dark + 0.005 * logphi0             
	e21      = e21dark + 0.004 * logphi0             

	S0       = 0.5 * (1. + tanh(120.*(100. * Irradiance - 0.1))) 
	Ephoton  = 1.E9 * hc / wavelength          
	flux     = 1000. * Irradiance / Ephoton    
	F        = flux  * sigma_retinal / (wloss * 1000.) 
}
UNITSON
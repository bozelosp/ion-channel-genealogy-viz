NEURON {
	SUFFIX canmda 		
		
	
	USEION ca WRITE ica
	NONSPECIFIC_CURRENT i 
	RANGE g, i, mg, inmda, gnmda, iampa, gampa, itotal, irtype, Pca, P, f
	
	GLOBAL Area			
		
	
	EXTERNAL i2_nmda, g2_nmda, i2_ampa, g2_ampa, irtype_car
		
	}

UNITS {
	(mA) = (milliamp)
	(mV) = (millivolt)
	(uS) = (microsiemens)
}

PARAMETER {                     
        dt			(ms)
       	
		mg   = 1	(mM)        
			
		Area = 1.11e-8  (cm2)	
		k = 1e-06   (mA/nA)		

		P           (cm/s/uS) 	
}  


ASSIGNED {	
	ica (mA/cm2)	
	v (mV)          
	i (mA/cm2)		
	g (uS)          
	Pca (cm/s)		
	
	inmda (mA/cm2)	
	gnmda	(uS)	
	iampa	(mA/cm2)
	gampa	(uS)	
	itotal (mA/cm2)	
	irtype (mA/cm2)	
	f               
}

INITIAL {

	P  = (1-exp(-65*-0.0755))/(10*Area*14564*(50e-09-(2e-03*exp(-65*-0.0755))))*k	
}


BREAKPOINT {
	g = g2_nmda	
	Pca = P*g	
	ica = Pca*14564*v*(50e-09-(2e-03*exp(v*-0.0755)))/(1-exp(v*-0.0755))*1/(1+(exp(0.08(/mV) * -v)*(mg / 0.69)))	
	i = -Pca*14564*v*(50e-09-(2e-03*exp(v*-0.0755)))/(1-exp(v*-0.0755))*1/(1+(exp(0.08(/mV) * -v)*(mg / 0.69)))	
	
	
	
	

	gnmda=g2_nmda*1/(1+(exp(0.08(/mV) * -v)*(mg / 0.69)))	
	gampa=g2_ampa	
	inmda=-i2_nmda	
	iampa=-i2_ampa	
	
	irtype=irtype_car	
	itotal=i2_nmda+i2_ampa+irtype_car		
	f=i/inmda		
	
}
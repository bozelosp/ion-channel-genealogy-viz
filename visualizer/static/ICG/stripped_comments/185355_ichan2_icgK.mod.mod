UNITS {
        (mA) =(milliamp)
        (mV) =(millivolt)
        (uF) = (microfarad)
	(molar) = (1/liter)
	(nA) = (nanoamp)
	(mM) = (millimolar)
	(um) = (micron)
	FARADAY = 96520 (coul)
	R = 8.3134	(joule/degC)
}
 
NEURON { 

SUFFIX ichan2



USEION k READ ek WRITE ik  				




RANGE gnatbar, gkfbar, gksbar				
RANGE gl, el, ina, ik, il, ggabaa, igabaa, egabaa	
GLOBAL ek
}



 


PARAMETER {						
        
        
        

        gnatbar (mho/cm2)   				
        ena  	(mV)	
		
	gkfbar =1.0	(mho/cm2)				
	gksbar = 1.0 (mho/cm2)	                        
        ek     	(mV)                      

	gl 	(mho/cm2)    				
 	el 	(mV)

	ggabaa 	(mho/cm2)    				
 	egabaa 	(mV)
}



STATE {
	m h nf ns
}




ASSIGNED {		

        v (mV) 
        celsius (degC)
        dt (ms) 
	

        gna (mho/cm2) 					
        ina (mA/cm2)
	
        gkf (mho/cm2)					
        gks (mho/cm2)
	ik (mA/cm2)

	il (mA/cm2)					

	igabaa (mA/cm2)					

	minf hinf nfinf nsinf				
 	mtau (ms) htau (ms) nftau (ms) nstau (ms)	
	mexp hexp nfexp nsexp
	q10
} 


BREAKPOINT {
	SOLVE states					
        gna = gnatbar*m*m*m*h  			
        gkf = gkfbar*nf*nf*nf*nf
        gks = gksbar*ns*ns*ns*ns

        
       	ik = gkf*(v-ek) + gks*(v-ek)
	
	
}
 

 

INITIAL {
	trates(v)
	
	m = minf
	h = hinf
        nf = nfinf
        ns = nsinf
	
	VERBATIM
	return 0;
	ENDVERBATIM
}


PROCEDURE states() {	
        trates(v)	
        m = m + mexp*(minf-m)
        h = h + hexp*(hinf-h)
        nf = nf + nfexp*(nfinf-nf)
        ns = ns + nsexp*(nsinf-ns)
        VERBATIM
        return 0;
        ENDVERBATIM
}





PROCEDURE rates(v) {  
        LOCAL  alpha, beta, sum
        q10 = 3^((celsius - 6.3)/10)
                
	alpha = -0.3*vtrap((v+60-17),-5)		
	beta = 0.3*vtrap((v+60-45),5)			
	sum = alpha+beta        
	mtau = 1/sum      minf = alpha/sum
                
	alpha = 0.23/exp((v+60+5)/20)			
	beta = 3.33/(1+exp((v+60-47.5)/-10))		
	sum = alpha+beta
	htau = 1/sum 
        hinf = alpha/sum 


             
        alpha = -0.028*vtrap((v+65-35),-6)		
	beta = 0.1056/exp((v+65-10)/40)			
	sum = alpha+beta        			
	nstau = 1/sum      nsinf = alpha/sum		
            
        alpha = -0.07*vtrap((v+65-47),-6)		
	beta = 0.264/exp((v+65-22)/40)			
	sum = alpha+beta        
	nftau = 1/sum      nfinf = alpha/sum
}


PROCEDURE trates(v) {  
	LOCAL tinc
        
	
                           
	rates(v)	
			
			

        tinc = -dt * q10
        mexp = 1 - exp(tinc/mtau)
        hexp = 1 - exp(tinc/htau)
	nfexp = 1 - exp(tinc/nftau)
	nsexp = 1 - exp(tinc/nstau)
}


FUNCTION vtrap(x,y) {  
        if (fabs(x/y) < 1e-6) {
                vtrap = y*(1 - x/y/2)
        }else{  
                vtrap = x/(exp(x/y) - 1)
        }
}
UNITS {
        (mA) =		(milliamp)
        (mV) =		(millivolt)
        (uF) = 		(microfarad)
	(molar) = 	(1/liter)
	(nA) = 		(nanoamp)
	(mM) = 		(millimolar)
	(um) = 		(micron)
	FARADAY = 96520 (coul)
	R = 8.3134	(joule/degC)
}
 
NEURON { 
SUFFIX lca
USEION lca READ elca WRITE ilca VALENCE 2 
RANGE  glca
RANGE glcabar
RANGE einf, etau, ilca
}
 
INDEPENDENT {t FROM 0 TO 100 WITH 100 (ms)}
 
PARAMETER {
        v (mV) 
        celsius = 6.3 (degC)
        dt (ms) 
	glcabar (mho/cm2)
}
 
STATE {
	e
}
 
ASSIGNED {
	glca (mho/cm2)
	ilca (mA/cm2)
	elca (mV)

	einf 
	etau (ms) 
	eexp      
} 

BREAKPOINT {
	SOLVE states
        glca = glcabar*e*e
	ilca = glca*(v-elca)
}
 
UNITSOFF
 
INITIAL {
	trates(v)
	e = einf
}

PROCEDURE states() {	
        trates(v)	
	e = e + eexp*(einf-e)
        VERBATIM
        return 0;
        ENDVERBATIM
}
 
LOCAL q10

PROCEDURE rates(v) {  
                      
        LOCAL  alpha, beta, sum
        q10 = 3^((celsius - 6.3)/10) 
                
        alpha = -15.69*vtrap(v-81.5,-10)	
	beta = 0.29*exp(-v/10.86)	
	sum = alpha+beta        
	etau = 1/sum      einf = alpha/sum
                
}
 
PROCEDURE trates(v) {  
                      
	LOCAL tinc
        TABLE  einf, eexp, etau
	DEPEND dt, celsius FROM -100 TO 100 WITH 200
                           
	rates(v)	
		
		

	tinc = -dt * q10
	eexp = 1 - exp(tinc/etau)
}
 
FUNCTION vtrap(x,y) {  
        if (fabs(x/y) < 1e-6) {
                vtrap = y*(1 - x/y/2)
        }else{  
                vtrap = x/(exp(x/y) - 1)
        }
}
 
UNITSON
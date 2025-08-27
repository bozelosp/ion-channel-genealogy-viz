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
 
? interface 
NEURON { 
SUFFIX ichan 



USEION na READ ena WRITE ina
USEION k READ ek WRITE ik
RANGE  gnat, gkf
RANGE gnatbar, gkfbar
RANGE gl, el
RANGE minf, mtau, hinf, htau, nfinf, nftau
}
 
INDEPENDENT {t FROM 0 TO 100 WITH 100 (ms)}
 
PARAMETER {
        v (mV) 
        celsius = 6.3 (degC)
        dt (ms) 
        ena  (mV)
	gnatbar = 0.0 (mho/cm2)   
        ek  (mV)
	gkfbar = 1.0 (mho/cm2)
	gl = 0.0(mho/cm2)    
 	el (mV)
}
 
STATE {
	m h nf
}
 
ASSIGNED {
         
        gnat (mho/cm2) 
        gkf (mho/cm2)


        ina (mA/cm2)
        ik (mA/cm2)

	il (mA/cm2)

	minf hinf nfinf
 	mtau (ms) htau (ms) nftau (ms)
	mexp hexp nfexp
} 

? currents
BREAKPOINT {
	SOLVE states
        gnat = gnatbar*m*m*m*h  
        ina = gnat*(v - ena)
        gkf = gkfbar*nf*nf*nf*nf
        ik = gkf*(v-ek)
	il = gl*(v-el)
}
 
UNITSOFF
 
INITIAL {
	trates(v)
	
	m = minf
	h = hinf
      nf = nfinf
	
	VERBATIM
	return 0;
	ENDVERBATIM
}

? states
PROCEDURE states() {	
        trates(v)	
        m = m + mexp*(minf-m)
        h = h + hexp*(hinf-h)
        nf = nf + nfexp*(nfinf-nf)
        VERBATIM
        return 0;
        ENDVERBATIM
}
 
LOCAL q10

? rates
PROCEDURE rates(v) {  
                      
        LOCAL  alpha, beta, sum
       q10 = 3^((celsius - 6.3)/10)
                
	alpha = -0.3*vtrap((v+65-25),-5)
	beta = 0.3*vtrap((v+65-53),5)
	sum = alpha+beta        
	mtau = 1/sum      minf = alpha/sum
                
	alpha = 0.23/exp((v+65-3)/20)
	beta = 3.33/(1+exp((v+65-55.5)/-10))
	sum = alpha+beta
	htau = 1/sum 
        hinf = alpha/sum 
             
        alpha = -0.07*vtrap((v+65-47),-6)
	beta = 0.264/exp((v+65-22)/40)
	sum = alpha+beta        
	nftau = 1/sum      nfinf = alpha/sum
	
}
 
PROCEDURE trates(v) {  
                      
	LOCAL tinc
        TABLE minf, mexp, hinf, hexp, nfinf, nfexp, mtau, htau, nftau
	DEPEND dt, celsius FROM -100 TO 100 WITH 200
                           
	rates(v)	
		
		

	       tinc = -dt * q10
        mexp = 1 - exp(tinc/mtau)
        hexp = 1 - exp(tinc/htau)
	nfexp = 1 - exp(tinc/nftau)
}
 
FUNCTION vtrap(x,y) {  
        if (fabs(x/y) < 1e-6) {
                vtrap = y*(1 - x/y/2)
        }else{  
                vtrap = x/(exp(x/y) - 1)
        }
}
 
UNITSON
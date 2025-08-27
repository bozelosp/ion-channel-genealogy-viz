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
SUFFIX ichan2 


USEION k READ ek WRITE ik

RANGE  gkf, gks
RANGE gkfbar, gksbar
RANGE gl, el
RANGE nfinf, nftau, ikf, nsinf, nstau, iks
}
 
INDEPENDENT {t FROM 0 TO 100 WITH 100 (ms)}
 
PARAMETER {
        v (mV) 
        celsius = 6.3 (degC)
        dt (ms) 
        ekf  (mV)
	gkfbar = 1.0 (mho/cm2)
        eks  (mV)
	gksbar = 1.0 (mho/cm2)
	gl (mho/cm2)    
 	el (mV)
}
 
STATE {
	nf ns
}
 
ASSIGNED {
	ik (mA/cm2)
        ek (mV) 
        gkf (mho/cm2)
        gks (mho/cm2)

        ikf (mA/cm2)
        iks (mA/cm2)


	il (mA/cm2)

	nfinf nsinf
 	nftau (ms) nstau (ms)
	nfexp nsexp
} 

? currents
BREAKPOINT {
	SOLVE states
        gkf = gkfbar*nf*nf*nf*nf
        ikf = gkf*(v-ekf)
        gks = gksbar*ns*ns*ns*ns
        iks = gks*(v-eks)
	ik = ikf

	il = gl*(v-el)
}
 
UNITSOFF
 
INITIAL {
	trates(v)
	
      nf = nfinf
      ns = nsinf

}

? states
PROCEDURE states() {	
        trates(v)	
        nf = nf + nfexp*(nfinf-nf)
        ns = ns + nsexp*(nsinf-ns)
}
 
LOCAL q10

? rates
PROCEDURE rates(v) {  
                      
        LOCAL  alpha, beta, sum
       q10 = 3^((celsius - 6.3)/10)
                
             
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
        TABLE nfinf, nfexp, nsinf, nsexp, nftau, nstau
	DEPEND dt, celsius FROM -100 TO 100 WITH 200
                           
	rates(v)	
		
		

	       tinc = -dt * q10
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
 
UNITSON
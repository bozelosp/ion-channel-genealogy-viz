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
SUFFIX hyperde1 




NONSPECIFIC_CURRENT i
RANGE  ghyf, ghys, ghyhtf, ghyhts
RANGE ghyfbar, ghysbar, ghyhtfbar, ghyhtsbar
RANGE hyfinf, hysinf, hyftau, hystau
RANGE hyhtfinf, hyhtsinf, hyhtftau, hyhtstau
GLOBAL eh
}
 
INDEPENDENT {t FROM 0 TO 100 WITH 100 (ms)}
 
PARAMETER {
      v (mV) 
      celsius = 6.3 (degC)
      dt (ms) 

	ghyfbar = 1.0 (mho/cm2)
	ghysbar = 1.0 (mho/cm2)
	ehyf (mV)
	ehys (mV)
	ghyhtfbar (mho/cm2)
	ghyhtsbar (mho/cm2)
	ehyhtf (mV)
	ehyhts (mV)
}
 
STATE {
	hyf hys hyhtf hyhts
}
 
ASSIGNED {
         
  
	ghyf (mho/cm2)
 	ghys (mho/cm2)

	ghyhtf (mho/cm2)
	ghyhts (mho/cm2)

	i (mA/cm2)
	eh (mV)
  
	ihyf (mA/cm2)
	ihys (mA/cm2)
	ihyhtf (mA/cm2)
	ihyhts (mA/cm2)

	hyfinf hysinf hyhtfinf hyhtsinf
 	hyftau (ms) hystau (ms) hyhtftau (ms) hyhtstau (ms)
	hyfexp hysexp hyhtfexp hyhtsexp     
} 

? currents
BREAKPOINT {

	SOLVE states

	
	
	ghys = ghysbar * hys*hys
	ihys = ghys * (v-eh)
	i = ihys

	
	
	
	
		
		}
 
UNITSOFF
 
INITIAL {
	trates(v)
	
	hyf = hyfinf
      hys = hysinf
	hyhtf = hyhtfinf
	hyhts = hyhtsinf
	VERBATIM
	return 0;
	ENDVERBATIM
}

? states
PROCEDURE states() {	
        trates(v)	
        
        hyf = hyf + hyfexp*(hyfinf-hyf)
        hys = hys + hysexp*(hysinf-hys)
	  hyhtf = hyhtf + hyhtfexp*(hyhtfinf-hyhtf)
	  hyhts = hyhts + hyhtsexp*(hyhtsinf-hyhts)

        VERBATIM
        return 0;
        ENDVERBATIM
}
 
LOCAL q10

? rates
PROCEDURE rates(v) {  
                      
        LOCAL  alpha, beta, sum
       q10 = 3^((celsius - 6.3)/10)
       
	
	hyfinf =  1 / (1 + exp( (v+81+3.33)/10 ))
	hyftau = 14.9 + 14.1 / (1+exp(-(v+95.2)/0.5))

	
	hysinf =  1 / (1 + exp( (v+81+3.33)/10 ))
	hystau = 80 + 172.7 / (1+exp(-(v+59.3)/-0.83))

		
	hyhtfinf =  1 / (1 + exp( (v+77+3.33)/10 ))
	hyhtftau = 23.2 + 16.1 / (1+exp(-(v+91.2)/0.83))

		
	hyhtsinf =  1 / (1 + exp( (v+77+3.33)/10 ))
	hyhtstau = 227.3 + 170.7*exp(-0.5*((v+80.4)/11)^2)
}
 
PROCEDURE trates(v) {  
                      
	LOCAL tinc
      TABLE hyfinf, hyhtfinf, hyfexp, hyhtfexp, hyftau, hyhtftau, 
		hysinf, hyhtsinf, hysexp, hyhtsexp, hystau, hyhtstau	
	DEPEND dt, celsius FROM -120 TO 100 WITH 220
                           
	rates(v)	
		
		

	       tinc = -dt * q10
        
        hyfexp = 1 - exp(tinc/hyftau)
	  hysexp = 1 - exp(tinc/hystau)
	  hyhtfexp = 1 - exp(tinc/hyhtftau)
	  hyhtsexp = 1 - exp(tinc/hyhtstau)
}
 
FUNCTION vtrap(x,y) {  
        if (fabs(x/y) < 1e-6) {
                vtrap = y*(1 - x/y/2)
        }else{  
                vtrap = x/(exp(x/y) - 1)
        }
}
 
UNITSON
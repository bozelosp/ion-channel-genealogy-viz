NEURON {
	SUFFIX calH
	USEION ca READ eca WRITE ica
        RANGE gcalbar, m, h,sh
	RANGE inf, fac, tau, mytau, myhtau
}

UNITS {
	(mA) = (milliamp)
	(mV) = (millivolt)
}



PARAMETER {          
        v               (mV)
        sh=0              (mV)

        celsius = 34	(degC)
	dt              (ms)
        gcalbar = 0     (mho/cm2) 
	eca = 140       (mV)      
      mytau=3.6 (ms)
      myhtau=29 (ms)  

        }

STATE {	m h }                     

ASSIGNED {
	ica (mA/cm2)
      inf[2]
	tau[2]

        
}


INITIAL {
      m = 0    
	h = 1    
	rate(v)
	
}

BREAKPOINT {
	SOLVE state METHOD cnexp
	ica = gcalbar*m*m*m*h*(v - eca)

}



DERIVATIVE state {  
        rate(v)
        m' = (inf[0]-m)/tau[0]
	  h' = (inf[1]-h)/tau[1]

}

PROCEDURE rate(v (mV)) { 
       FROM i=0 TO 1 {
		tau[i] = vartau(v,i)
		inf[i] = varss(v,i)
	}

     
	
}


FUNCTION varss(v, i) {
	if (i==0) { 
             varss = 1 / (1 + exp((v+37-sh)/(-1)))  
	}
	else if (i==1) { 
             varss = 1 / (1 + exp((v+41-sh)/(0.5))) 
	}
}

FUNCTION vartau(v, i) {
	if (i==0) {
          vartau = mytau  
         

        }
	else if (i==1) {
           vartau = myhtau   
  
        }
}
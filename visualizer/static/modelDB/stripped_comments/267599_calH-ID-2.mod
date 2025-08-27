NEURON {
	SUFFIX calH
	USEION ca READ eca WRITE ica
        RANGE gcal,gcalbar, m, h,ica
	RANGE inf, fac, tau
}

UNITS {
	(mA) = (milliamp)
	(mV) = (millivolt)
}



PARAMETER {          
        v               (mV)
        celsius = 34	(degC)
	dt              (ms)
        gcalbar = 0     (mho/cm2) 
        gcal            (mho/cm2) 
	eca = 140       (mV)      
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
	gcal = gcalbar*m*m*m*h

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
	LOCAL Arg
    if (i==0) {
        Arg=(v+37)/(-1)
        
        if (Arg<-50) {varss=1}              
        else if (Arg>50) {varss=0}
        else {varss=1/(1+exp(Arg))}
    } else if (i==1) {
        Arg=(v+41)/(0.5)
        
        if (Arg<-50) {varss=1}              
        else if (Arg>50) {varss=0}
        else {varss=1/(1+exp(Arg))}
    }
}

FUNCTION vartau(v, i) {
	if (i==0) {
          vartau = 3.6  
         

        }
	else if (i==1) {

           vartau = 29   
        }
}
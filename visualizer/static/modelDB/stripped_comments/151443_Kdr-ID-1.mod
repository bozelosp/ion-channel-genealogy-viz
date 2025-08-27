UNITS {
        (mA) = (milliamp)
        (mV) = (millivolt)
		(S) = (siemens)
}
 
NEURON {
        SUFFIX Kdr
        USEION k READ ek WRITE ik
        RANGE gkmax, gk
        RANGE ninf, ntau
}
 
PARAMETER { 
        gkmax = .5 (S/cm2)	  <0,1e9>
}
 
STATE {
        n
}
 
ASSIGNED {
        v (mV)
        ek (mV)

		gk (S/cm2)
        ik (mA/cm2)
        ninf
		ntau (ms)
}
 
BREAKPOINT {
        SOLVE states METHOD cnexp
        gk = gkmax*n*n*n*n
		ik = gk*(v - ek)
} 
 
INITIAL {
	rates(v)
	n = ninf
}

DERIVATIVE states {
        rates(v)
        n' = (ninf-n)/ntau
}
 
PROCEDURE rates(v(mV)) {  
		
        LOCAL  alpha, beta, sum

UNITSOFF
                
        ntau = 5/(exp((50+v)/40)+exp(-(50+v)/50))
        ninf = 1/(1+exp((38+v)/-15))
}
 
FUNCTION vtrap(x,y) {  
        if (fabs(x/y) < 1e-6) {
                vtrap = y*(1 - x/y/2)
        }else{
                vtrap = x/(exp(x/y) - 1)
        }
}
 
UNITSON
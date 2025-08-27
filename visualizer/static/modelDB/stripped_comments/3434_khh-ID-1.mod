UNITS {
        (S)  = (siemens)
        (mV) = (millivolt)
        (mA) = (milliamp)
}
 
NEURON {
        SUFFIX khh
        USEION k READ ek WRITE ik
        RANGE  gkbar, gk, ik, alphaVHalf, alphaK, alpha0, betaVHalf, betaK, beta0, q10, T0
}


 
PARAMETER {
        gkbar = .036 (S/cm2)	<0,1e9>
        alpha0 = 0.1 (/ms)	    <0,1e9>
        alphaVHalf = -55 (mV)
        alphaK = 10.0 (mV)		<0,1e9>
        beta0 = 0.125 (/ms)	    <0,1e9>
        betaVHalf = -65 (mV)
        betaK = 80.0 (mV)	    <0,1e9>
        q10	= 3.0       	    <0,1e9>
        T0	= 6.3 (degC)
}
 
ASSIGNED {
        gk (mho/cm2)
        ik (mA/cm2)
        ratefac
        alphaKInv (/mV)
        betaKInv (/mV)
        v (mV)
        celsius (degC)
        ek (mV)
}

STATE { n }
 
BREAKPOINT {
        SOLVE states METHOD cnexp
        gk = gkbar*n*n*n*n
        ik = gk*(v - ek)      
}
 
INITIAL {
        n = alpha(v)/(alpha(v) + beta(v)) 
        ratefac = q10^((celsius - T0)/10(degC)) 
        alphaKInv = 1/alphaK 
        betaKInv = 1/betaK
}

DERIVATIVE states {
        n' = ((1-n)*alpha(v) - n*beta(v))*ratefac
}
 
FUNCTION alpha(Vm (mV)) (/ms) { 
         LOCAL x
         x = (Vm-alphaVHalf)*alphaKInv
         if (fabs(x) > 1e-6) {            
                alpha = alpha0*x/(1-exp(-x))
        }else{
                alpha = alpha0/(1 - 0.5*x)
        }
}
 
FUNCTION beta(Vm (mV)) (/ms) { 
         beta = beta0*exp((-Vm+betaVHalf)*betaKInv)
}
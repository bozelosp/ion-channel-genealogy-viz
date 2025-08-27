UNITS {
        (S)  = (siemens)
        (mV) = (millivolt)
        (mA) = (milliamp)
}
 
NEURON {
        SUFFIX khva
        USEION k READ ek WRITE ik
        RANGE  gkbar, gk, ik, alphaVHalf, alphaK, alpha0, betaVHalf, betaK, beta0, q10, T0
}


 
PARAMETER {
        gkbar = 0.04 (S/cm2)	<0,1e9>
        alpha0 = 0.11 (/ms)		<0,1e9>
        alphaVHalf = -19 (mV)
        alphaK = 9.10 (mV)		<0,1e9>
        beta0 = 0.103 (/ms)		<0,1e9>
        betaVHalf = -19 (mV)
        betaK = 20.0 (mV)		<0,1e9>
        q10	= 2.0				<0,1e9>
        T0	= 23 (degC)			<0,1e9>
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
        gk = gkbar*n 
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
       alpha = alpha0*exp((Vm-alphaVHalf)*alphaKInv) 
}
 
FUNCTION beta(Vm (mV)) (/ms) { 
        beta  = beta0*exp((-Vm+betaVHalf)*betaKInv)
}
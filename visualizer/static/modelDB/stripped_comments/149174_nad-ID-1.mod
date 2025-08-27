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
SUFFIX nad 
USEION nat READ enat WRITE inat VALENCE 1
RANGE gnat
RANGE gnatbar
RANGE inat
RANGE alphaD, betaD
}
 
INDEPENDENT {t FROM 0 TO 100 WITH 100 (ms)}
 
PARAMETER {
    v (mV) 
    celsius = 6.3 (degC)
    dt (ms) 
    enat  (mV)
	gnatbar (mho/cm2)
    alphaD (/ms)
    betaD (/ms)
}
 
STATE {
O C1 C2 C3 I I1 I2 I3 ID I1D I2D I3D
}
 
KINETIC scheme1 {
rates(v)

~ O  <-> C1  (3*betam, 1*alpham)
~ O  <-> I   (1*betah, 1*alphah)
~ C1 <-> C2  (2*betam, 2*alpham)
~ C2 <-> C3  (1*betam, 3*alpham)
~ I  <-> I1  (3*betam, 1*alpham)
~ I1 <-> I2  (2*betam, 2*alpham)
~ I2 <-> I3  (1*betam, 3*alpham)
~ C1 <-> I1  (1*betah, 1*alphah)
~ C2 <-> I2  (1*betah, 1*alphah)
~ C3 <-> I3  (1*betah, 1*alphah)
~ I  <-> ID  (betaD,   alphaD)
~ I1 <-> I1D (betaD,   alphaD)
~ I2 <-> I2D (betaD,   alphaD)
~ I3 <-> I3D (betaD,   alphaD)

CONSERVE O+C1+C2+C3+I+I1+I2+I3+ID+I1D+I2D+I3D = 1
}


ASSIGNED {
        gnat (mho/cm2) 
        inat (mA/cm2)
        alpham (/ms)
        alphah (/ms)
        betam (/ms)
        betah (/ms)
} 

? currents
BREAKPOINT {
	SOLVE scheme1 METHOD sparse
        gnat = gnatbar*O  
        inat = gnat*(v - enat)
}
 
UNITSOFF
 
INITIAL {
	rates(v)
	SOLVE scheme1 STEADYSTATE sparse
}

LOCAL q10

? rates
PROCEDURE rates(v) {  
                      
       q10 = 3^((celsius - 6.3)/10)
                
	alpham = -0.3*vtrap((v+60-17),-5)
	betam = 0.3*vtrap((v+60-45),5)
                
	alphah = 0.23/exp((v+60+5)/20)
	betah = 3.33/(1+exp((v+60-47.5)/-10))
}
 
FUNCTION vtrap(x,y) {  
        if (fabs(x/y) < 1e-6) {
                vtrap = y*(1 - x/y/2)
        }else{  
                vtrap = x/(exp(x/y) - 1)
        }
}
 
UNITSON
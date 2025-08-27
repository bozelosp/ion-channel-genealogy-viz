UNITS {
       (molar) = (1/liter)
       (S)  = (siemens)
       (mA) = (milliamp)
       (mV) = (millivolt)
       (mM) = (millimolar)
        F = (faraday)  (coulomb)
        R = (mole k)   (mV-coulomb/degC)
       
}
 
NEURON {
        SUFFIX cachan
        USEION ca READ cai WRITE ica
        RANGE  gcalbar,gcanbar,gcahvabar,ica,ical,icahva,ican,kml,kmn, mhalf, mslope
        GLOBAL dlinf
}
 
INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}
 
PARAMETER {
        v (mV)
        dt (ms)
        cai   (mM)
        celsius =  35.0      (degC)
        gcahvabar =  0.0e-6 (S/cm2)
        gcalbar =  11.196e-6  (S/cm2)
        kmn = 0.0001   (mM)
        kml = 0.00045  (mM)
        mhalf=-35
        mslope = 7
        eca = 120 (mV)
        cao = 2.0 (mM)
        
}
 
STATE {
        dl
}
 
ASSIGNED {
        ica (mA/cm2)
        ical (mA/cm2)
        dlinf
 }
 
BREAKPOINT {
        SOLVE states METHOD cnexp
        ical = gcalbar*dl*(v - eca)
        
        ica  = ical
}
 
UNITSOFF
 
INITIAL {
        dl = boltz(v,mhalf,mslope)
        
        
}

DERIVATIVE states {  
LOCAL dlinf,dhvainf,fhvainf,dltau,dhvatau,fhvatau
        dlinf = boltz(v,mhalf,mslope)
        
        
        dltau = gaussian(v,9.0,25.0,70.0,0.30)
        
        
        dl'  = (dlinf-dl)/dltau
        
        
}
 
 
FUNCTION gaussian(v,a,b,c,d) {
        LOCAL arg
        arg= a*exp(-(c+v)*(v+c)/(b*b)) +d
        gaussian = arg
}
 
 
FUNCTION boltz(x,y,z) {
               LOCAL arg
                arg= -(x-y)/z
                if (arg > 50) {boltz = 0}
                else {if (arg < -50) {boltz = 1}
                else {boltz = 1.0/(1.0 + exp(arg))}}
}

 
UNITSON
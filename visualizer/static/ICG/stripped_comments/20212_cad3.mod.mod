INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}

NEURON {
        SUFFIX cad3
        USEION ca READ ica, cai WRITE cai
        RANGE m
        GLOBAL depth,kt,kd,cainf,taur
}

UNITS {
        (molar) = (1/liter)                     
        (mM)    = (millimolar)
        (um)    = (micron)
        (mA)    = (milliamp)
        (msM)   = (ms mM)
        FARADAY = (faraday) (coulomb)
}






PARAMETER {
        depth   = .1    (um)            
        
        taur    = 200   (ms)            
        cainf   = 1e-8  (mM)
        cai		(mM)
        kt      = 1     (mM/ms)         
       
        kd      = 5e-4  (mM)            
}

STATE {
        m             (mM) 
}

INITIAL {
        m =  cainf 
}

ASSIGNED {
        ica             (mA/cm2)
        drive_channel   (mM/ms)
        drive_pump      (mM/ms)
}
        
BREAKPOINT {
        SOLVE state METHOD derivimplicit
}

DERIVATIVE state { 

        drive_channel =  - (10000) * ica / (2 * FARADAY * depth)

        if (drive_channel <= 0. ) { drive_channel = 0. } 



        drive_pump = -kt * m / (m + kd )            

         m' = drive_channel + drive_pump + (cainf-m)/taur
         cai = m
         
}
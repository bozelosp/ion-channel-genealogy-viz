NEURON {
        SUFFIX cadyn
        USEION ca READ cai,ica WRITE cai 
        RANGE CAF, tca, cai
}

UNITS {
        (mM) = (milli/liter)
        (mA) = (milliamp)
        F    = (faraday) (coul)
}

PARAMETER {
        tca= 70 (ms)           
        cainf= 50e-6   (mM)      
        dep= 2e-4 (micron)     

}

ASSIGNED {
        ica     (mA/cm2)
        diam    (micron)
        A       (/coul/cm)
        CAF     ()
}

STATE { cai (mM) }

BREAKPOINT { 
        SOLVE states METHOD cnexp
}

DERIVATIVE states {
         cai'= -A*CAF*ica - (cai-cainf)/tca
}

INITIAL {
        A =(1e4)/(F*dep)
}
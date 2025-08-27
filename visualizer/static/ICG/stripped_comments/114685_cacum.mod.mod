NEURON {
    SUFFIX cacum
    USEION ca READ ica WRITE cai
    GLOBAL depth, tau, cai0
}

UNITS {
    (mM) = (milli/liter)
    (mA) = (milliamp)
    F = (faraday) (coulombs)	
}

PARAMETER {
    depth = 200 (nm)      
    tau = 10 (ms)
    cai0 = 5e-5 (mM)
    cai0_ca_ion
}

ASSIGNED {
    ica (mA/cm2)
}

STATE {
    cai (mM)
}

INITIAL {
    cai=cai0
}

BREAKPOINT {
    SOLVE integrate METHOD derivimplicit
}

DERIVATIVE integrate {
    cai' = -ica/depth/F/2 * (1e7) + (cai0 - cai)/tau
}
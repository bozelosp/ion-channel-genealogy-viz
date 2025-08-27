UNITS {
   (mv) = (millivolt)
   (mA) = (milliamp)
}

NEURON {
    SUFFIX NaL
    USEION na READ ena,nai WRITE ina
    RANGE gna,inaL
    GLOBAL gmax_naL
}

PARAMETER {
    v (mV)
    gna = 2.8e-5   (mho/cm2)
    inaL = 0.0      (mA/cm2)
    ena
    nai
    celsius
}

ASSIGNED { 
    ina (mA/cm2)
    gmax_naL
}

BREAKPOINT {
    ina	= gna*gmax_naL*(v-ena)
    inaL = ina
}

INITIAL {
    gmax_naL = 1
}
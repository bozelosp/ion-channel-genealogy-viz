NEURON {
    SUFFIX myions
    USEION k WRITE ko, ki
    USEION na WRITE nao, nai
}

UNITS {
    (molar)	= (1/liter)
    (mM)	= (millimolar)
}

ASSIGNED {
    nai (mM)
    nao (mM)
    ki (mM)
    ko (mM)
}

BREAKPOINT {
    nai = nai
    nao = nao
    ki = ki
    ko = ko
}
INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}

NEURON {

        SUFFIX NOISE
        NONSPECIFIC_CURRENT i
        RANGE imax

}

ASSIGNED {

}

UNITS {
        (nA) = (nanoamp)
}

PARAMETER {
        imax=0          (umho)

}

INITIAL {


}

ASSIGNED { i (nA) }

BREAKPOINT {

SOLVE dum       

}

PROCEDURE dum() {

        i = (scop_random()-0.5)*imax

        
        
        VERBATIM
                return 0;
        ENDVERBATIM

}
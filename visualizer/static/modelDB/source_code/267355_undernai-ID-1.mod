TITLE undersampling of nai

COMMENT
	undersampling of the simulated nai to match the experimental speed 
ENDCOMMENT


NEURON {
        SUFFIX undernai     
	USEION na READ nai
	
    RANGE nau,b, naund,tr
 
}

PARAMETER {
tSta=0.1 
}


ASSIGNED {
 nai (milli/liter)
 nau (milli/liter)
 naund (milli/liter)
 b
 tr 
 plp

}


INITIAL {

	b = 1
 plp=0.1
	nau = 10
	naund = 0



}


BREAKPOINT {
VERBATIM

plp=b*tSta;
naund=naund+nai;
tr=t;


if ( (tr < plp+0.005) && (tr > plp-0.005) ) {
	nau=naund/10;
naund=0;
b=b+1;

}




ENDVERBATIM
}

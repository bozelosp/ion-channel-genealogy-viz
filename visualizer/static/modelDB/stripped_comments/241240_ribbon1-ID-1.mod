VERBATIM

double du_dev0();
double dexp_dev(double);
int igeom_dev(double);
ENDVERBATIM


NEURON {
  POINT_PROCESS 	ribbon1
  RANGE 		i, del, xp, yp, zp, dist, tau, sw
  NONSPECIFIC_CURRENT 	i
}

PARAMETER {
  del	= -1	(ms)	
  xp	= 0		
  yp 	= 0
  zp 	= 0
  dist	= 0		
  e 	= 0   	(mV)	
  tau	= 1.3	(ms)	
  atau	= 0	(ms)	
  sw	= 1000	(pS)	
  mp	= 1.0		
  idebug = 0
}

ASSIGNED {
  v (mV)
  i (nanoamp)
}

STATE {
  a (pS)
  g (pS)
}

UNITS {
 (mV) 	= (millivolt)
 (pS) 	= (picosiemens)
 PI	= (pi) (1)
}

INITIAL {
  if( del >= 0 ){
    net_send( del, 3 )
  }
  net_send( 1.0e-10, 1 )	
  g = 0
}

BREAKPOINT {
  SOLVE state METHOD sparse
  i = g*( v - e )*(1e-06)	
}

KINETIC state {
  ~ a <-> g ( 1/tau, 0 )
  ~ g  ->   ( 1/tau )
}

NET_RECEIVE( weight (pS)) {
  LOCAL wait
  UNITSOFF
  if( idebug ) {
    printf( "Net_receive afhcsyn
    if( flag == 0 ){ printf( " weight %g", weight ) }
    printf( "\n" )
  }
  if( flag == 0 ){	
    a = a + weight * exp(1) * igeom_dev(1/mp)
  }
  if( flag == 1 && atau > 0 ){	
    wait = dexp_dev( 1(ms)/atau )
    net_send( wait, 2 )
  }
  if( flag == 2 && atau > 0 ){	
    a = a + sw *exp(1)* igeom_dev(1/mp)
    wait = dexp_dev( 1/atau )
    net_send( wait, 2 )
  }
  if( flag == 3 ){	
    a = a + sw *exp(1) * igeom_dev(1/mp)
  }
  UNITSON
}
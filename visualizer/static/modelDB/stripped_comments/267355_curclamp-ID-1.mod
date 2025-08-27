NEURON {
    SUFFIX curclamp
    RANGE delay1,amplitude1, duration1, icin, delay2, duration2, amplitude2
    NONSPECIFIC_CURRENT  i
    RANGE i, e, ef, g
}

PARAMETER {

    delay1 = 1.1 (ms)
    duration1 = 6 (ms)
    amplitude1=7 (nanoamp/cm2)
    icin=0 (nanoamp)
        delay2 = 1.1 (ms)
    duration2 = 6 (ms)
    amplitude2=7 (nanoamp/cm2)
   

}
UNITS {
    (mV) = (millivolt)
    (mA) = (milliamp)
    (nA) = (nanoamp)
    (pA) = (picoamp)
    (S)  = (siemens)
}
ASSIGNED {
    i   (nanoamp/cm2)
    v   (millivolt) 

}



BREAKPOINT {
    at_time(delay1)
    if (t < delay1) {
          i = icin
      }
    if (t>delay1 && t<delay1+duration1) {
          i = amplitude1
      }
     if (t>delay1+duration1 && t<delay2) {
          i = icin
      }
          if (t>delay2 && t<delay2+duration2) {
          i = amplitude2
      }
  }
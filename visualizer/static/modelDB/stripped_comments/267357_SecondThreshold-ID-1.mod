NEURON {
  POINT_PROCESS SecNet
  RANGE thresh
}
PARAMETER { thresh = -60 } 
ASSIGNED { v }
INITIAL {
  net_send(0, 1) 
}
NET_RECEIVE(w) {
  if (flag == 1) {
    WATCH (v > thresh) 2 
  } else if (flag == 2) {
    net_event(t) 
  }
}
NEURON	{ 
  POINT_PROCESS SpikeGenerator2
  RANGE y
  RANGE start, end,on
  RANGE time
}

PARAMETER {
	start		= 100 (ms)	
	end		= 1e10 (ms)	
}

ASSIGNED {
	y
	 
	event (ms)
	 

	on

	time[200] (ms)	
	indice
}


INITIAL {	
	 
	indice=0
	on = 1

	y = -90
	event = start 
	event_time()
	while (on == 1 && event < 0) {
		event_time()
	}
	if (on == 1) {
		net_send(event, 1)
	}
}	

PROCEDURE event_time() {
	event=event+time[indice] 
	indice=indice+1
	if (event > end) {
		on = 0
	}
}

NET_RECEIVE (w) {
	if (flag == 1 && on == 1) {
		y = 20
		net_event(t)
		event_time()
		net_send(event - t, 1)
		net_send(.1, 2)
	}
	if (flag == 2) {
		y = -90
	}
}
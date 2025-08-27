NEURON {
	POINT_PROCESS Counter
	RANGE count 
}

PARAMETER {
	count = 0
}

NET_RECEIVE(change) {
    count = count + change
}
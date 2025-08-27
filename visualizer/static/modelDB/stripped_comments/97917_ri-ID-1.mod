NEURON { SUFFIX nothing }

VERBATIM
#ifdef NRN_MECHANISM_DATA_IS_SOA
#define get_nnode(sec) _nrn_mechanism_get_nnode(sec)
#define get_node(sec, idx) _nrn_mechanism_get_node(sec, idx)
#else
#define get_nnode(sec) sec->nnode
#define get_node(sec, idx) sec->pnode[idx]
const char* secname();
#endif
ENDVERBATIM

PROCEDURE scale_connection_coef(x, factor) {
VERBATIM {
	Section* sec;
	Node* nd;
#if defined(t)
	NrnThread* _nt = nrn_threads;
#endif
	sec = chk_access();
	if (_lx <= 0. || _lx > 1.) {
		hoc_execerror("out of range, must be 0 < x <= 1", (char*)0);
	}
	
	
	if (_lx == 1.) {
		nd = get_node(sec, get_nnode(sec) - 1);
	}else{
		nd = get_node(sec, (int) (_lx*(double)(get_nnode(sec) - 1)));
	}
	
	NODEA(nd) *= _lfactor;
	NODEB(nd) *= _lfactor;
}
ENDVERBATIM
}
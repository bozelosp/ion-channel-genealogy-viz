NEURON {
 	SUFFIX nothing
 }
 
 VERBATIM
 #ifdef NRN_MECHANISM_DATA_IS_SOA
 #define get_child(sec) _nrn_mechanism_get_child(sec)
 #define get_sibling(sec) _nrn_mechanism_get_sibling(sec)
 #else
 #define get_child(sec) sec->child
 #define get_sibling(sec) sec->sibling
 #endif
 static void subtree(Section* sec, Symbol* sym) {
 	nrn_pushsec(sec);	
 	hoc_run_stmt(sym);	
 	nrn_popsec();
 	for (Section* child = get_child(sec); child; child = get_sibling(child)) {
 		subtree(child, sym);
 	}
 }
 #ifndef NRN_VERSION_GTEQ_8_2_0
 Section* chk_access();
 Symbol* hoc_parse_stmt();
 #endif
 ENDVERBATIM
 
 PROCEDURE subtree_traverse_all() {
   VERBATIM
   {
 	Symlist* symlist = (Symlist*)0;
 	subtree(chk_access(), hoc_parse_stmt(gargstr(1), &symlist));
 	
 	hoc_free_list(&symlist);
   }
   ENDVERBATIM
 }
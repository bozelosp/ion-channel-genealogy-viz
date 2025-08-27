// Function to prepare data based on active filters
function prepare_data(network_data, filter_state) {

	// Filter nodes and obtain their IDs
	const nodes_data = network_data.nodes.filter(node => filter_nodes(filter_state, node));

	const nodes_id = new Set(nodes_data.map(node => node.id));  // Using a Set for faster look-up

	// Get the similarity score filter value, default to 0 if not available
	let similarity_score = filter_state['similarity score'].filter_value ;
	
	// Filter links based on node IDs and similarity score
	const links_data = network_data.links.filter(d => {
		const source_id = typeof d.source === 'object' ? d.source.id : d.source;
		const target_id = typeof d.target === 'object' ? d.target.id : d.target;
	
		//console.log(`Checking link with source ${source_id} and target ${target_id}`);
	
		const condition_met = nodes_id.has(source_id) && nodes_id.has(target_id) ;
	
		if (condition_met && d.weight > similarity_score ) {
			return true;
		} else {
			return false;
		}
	});	

	return [ nodes_data, links_data ] ;
	
}

// Function to filter nodes based on filter state
function filter_nodes(filter_state, node) {

	// Check filter criteria and return the status
	const num_copies = filter_state['# copies']?.filter_value || 0;
	const icg_entry = filter_state['ICG entry']?.filter_value || false;
	const all_classes = filter_state['All']?.filter_value || false;

	const ion_class = node.original_model.ion_class;  // Accessing ion class from the nested object

	const is_active_ion_channel = Object.keys(filter_state).some(tag => 
		filter_state[tag]?.filter_value && tag === ion_class
	);

	return (
		node.num_of_identicals >= num_copies &&
		(!icg_entry || node.original_model.ICG) &&
		(all_classes || is_active_ion_channel)
	);

}
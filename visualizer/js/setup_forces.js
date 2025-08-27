function setup_forces(links_data, nodes_data, graph_width, graph_height) {
	// Constants
	const LINK_DISTANCE = 12;
	const CHARGE_STRENGTH_MULTIPLIER = -0.5;
	const CHARGE_STRENGTH_CONSTANT = -10.575;

	// Link Force
	const link_force = setup_link_force(links_data, LINK_DISTANCE);

	// Charge Force
	const charge_force = setup_charge_force(nodes_data, CHARGE_STRENGTH_MULTIPLIER, CHARGE_STRENGTH_CONSTANT);

	// Center Force
	const center_force = setup_center_force(graph_width, graph_height);

	return [ link_force, charge_force, center_force ];
}

function setup_link_force(links_data, distance) {
	return d3.forceLink(links_data)
		.distance(() => distance)
		.id(d => d.id);  // Ensure this ID corresponds to node IDs
}

function setup_charge_force(nodes_data, multiplier, constant) {
	return d3.forceManyBody()
		.strength(d => {
			const nodeData = nodes_data.find(node => node.id === d.id);
			if (nodeData) {
				return (multiplier * Math.pow(nodeData.num_of_identicals, 1.125) + constant) + 
					   (multiplier * Math.pow(nodeData.num_of_identicals, 1.1275) + constant);
			} else {
				return 0;
			}
		});
}

function setup_center_force(graph_width, graph_height) {
	const graph_x_center = graph_width / 2;
	const graph_y_center = graph_height / 2;

	return d3.forceCenter(graph_x_center, graph_y_center);
}

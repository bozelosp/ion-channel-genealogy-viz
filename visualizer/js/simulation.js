// Constants

const ALPHA_RESTART = 0.875;
const ALPHA_TARGET_RESTART = 0.125;
const FORCE_STRENGTH = 0.1675;

//d3.forceX(d => node_id_to_location[d.id][0] * graph_width).strength(FORCE_STRENGTH) : 
		

async function setup_simulation(split_var, nodes_data, links_data, node_id_to_location, link, node, label, charge_force, link_force, graph_width, graph_height, existingSimulation) {
  
	const force_x = split_var ? 
	d3.forceX(d => node_id_to_location[d.id][0] * graph_width).strength(FORCE_STRENGTH) :
	d3.forceX(graph_width / 2);
	
	const force_y = split_var ? 
	d3.forceY(d => node_id_to_location[d.id][1] * graph_height).strength(FORCE_STRENGTH) :
	d3.forceY(graph_height / 2);

	existingSimulation?.stop();
	existingSimulation=null;
	// Create a new simulation
	existingSimulation = d3.forceSimulation()
		.nodes(nodes_data)
		.force("charge", charge_force)
		.force("links", link_force)
		.force("x", force_x)
		.force("y", force_y)
		.on("tick", () => tick_actions(link, node, label));

	return existingSimulation;

} 

function tick_actions(link, node, label) {
	link
		.attr("x1", function(d) { return d.source.x; })
		.attr("y1", function(d) { return d.source.y; })
		.attr("x2", function(d) { return d.target.x; })
		.attr("y2", function(d) { return d.target.y; });

	node
		.attr("cx", function(d) { return d.x; })
		.attr("cy", function(d) { return d.y; });

	label
		.attr("x", function(d) { return d.x; })
		.attr("y", function(d) { return d.y; });
}
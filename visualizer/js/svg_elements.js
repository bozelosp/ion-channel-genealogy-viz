function draw_or_update_elements(nodes_data, links_data, local_group) {
	// Handle Links
	const link_update = local_group.selectAll("line").data(links_data, d => d.id);
	link_update.exit().remove(); // Exit
	const link_enter = link_update.enter().append("line");
	link_enter.attr("class", "graph-link");
	const link = link_enter.merge(link_update); // Merge
  
	// Handle Nodes
	const node_update = local_group.selectAll("circle").data(nodes_data, d => d.id);
	node_update.exit().remove(); // Exit
	const node_enter = node_update.enter().append("circle");
	node_enter.attr("r", d => calculate_node_radius(d))
	  .attr("class", "graph-node");
	const node = node_enter.merge(node_update); // Merge
  
	// Handle Labels
	const label_update = local_group.selectAll(".node-label").data(nodes_data, d => d.id);
	label_update.exit().remove(); // Exit
	const label_enter = label_update.enter().append("text")
	  .attr("class", "node-label")
	  .attr("text-anchor", "middle") // center the text
	  .attr("dy", ".35em"); // adjust position to align it with the circle vertically
	label_enter.text(d => d.label || '');
	const label = label_enter.merge(label_update); // Merge
  
	// Order the elements
	local_group.selectAll(".graph-link").lower();
	local_group.selectAll(".graph-node").raise();
	local_group.selectAll(".node-label").raise();
  
	return [node, link, label];
  }
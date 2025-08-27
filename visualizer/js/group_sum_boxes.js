

// I MIGHT NEED THIS BELOW, I MIGHT NOT
function create_summary_boxes(nodes_grouped_by_filter) {

    const group_summary_container = document.getElementById('group-summary-container');

    // Remove all existing child nodes
    while (group_summary_container.firstChild) {
        group_summary_container.removeChild(group_summary_container.firstChild);
    }

    Object.entries(nodes_grouped_by_filter).forEach(([group_key, group_nodes], index) => {
        const individual_summary_box = document.createElement('div');
        individual_summary_box.className = 'group-summary-box';
        individual_summary_box.id = `group-summary-box-${index + 1}`;
        individual_summary_box.setAttribute('y', 286 + (index * 206));

        const group_title = document.createElement('span');
        group_title.className = 'group-summary-box-title';
        group_title.textContent = group_key;

        const node_count = document.createElement('span');
        node_count.className = 'group-summary-box-node-count';
        node_count.textContent = `${group_nodes.length} nodes`;

        individual_summary_box.appendChild(group_title);
        individual_summary_box.appendChild(document.createElement('br'));
        individual_summary_box.appendChild(node_count);

        group_summary_container.appendChild(individual_summary_box);
    });

}


// Update the UI to reflect the new groupings
function update_node_group_summary_boxes(groups, node_radius_map) {
	const node_groups_container = document.getElementById('node-groups-box');
	[...node_groups_container.getElementsByClassName('node-group')].forEach(elem => elem.remove());
	let current_y_position = STARTING_Y_POSITION;
	groups.forEach(([key, , nodes], index) => {
		const label = key.replace(',', ': ').replace(/,/g, ' • ').replace('Supermodel 1 • Supermodel 2', 'Supermodel 1 & 2');
		const sanitizedLabel = label.replace(/: | • | entry/g, '').replace('1 & 2', '1&2');
		const new_group_element = document.createElement('div');
		new_group_element.className = 'node-group';
		new_group_element.id = `node-group-${index + 1}`;
		new_group_element.innerHTML = `<span style="color: #ffc42c;">${label}</span><br><span style="color: #ffffff; margin-top: 9px; display: inline-block;">${nodes.length} nodes</span>`;
		new_group_element.style.y = current_y_position;
		current_y_position += Y_POSITION_STEP;
		node_groups_container.appendChild(new_group_element);
		const group_properties = key.split(',');
		group_properties.forEach(gp => {
			node_radius_map[`${key} -> ${gp}`] = SMALL_NODE_RADIUS;
		});
		nodes_data_array.push({ id: label, label: sanitizedLabel, num_of_identicals: 150 });
		node_radius_map[label] = LARGE_NODE_RADIUS;
	});
}

/* OLD CODE */
// Define the function
function update_grouped_nodes_in_ui_v1() {

	// Run the 'group_nodes_by_selected_filters' function and destructure its return value into two variables
	let [grouped_node_positions, label_data_list] = group_nodes_by_selected_filters(nodes_data_array, filter_state, node_radius_map);

	// Get the container element in which node groups will be displayed
	const node_groups_container = document.getElementById('node-groups-box');

	// Clear previous sum groups
	for (let index = 1; index < 35; index++) {
		const elem = document.getElementById(`node-group-${index}`);
		if (elem) elem.remove();
		else break;
	}

	// Initialize the y-position for placing new group elements
	let current_y_position = 286;

	// Initialize an object to store the group labels
	const group_labels = {};

	// Loop through each group's label data
	label_data_list.forEach((label_data, index) => {

		// Create the label for the current group by replacing and formatting substrings
		let current_group_label = label_data[0].replace(',', ': ').replace(/,/g, ' • ').replace('Supermodel 1 • Supermodel 2', 'Supermodel 1 & 2');

		// Create a new HTML div element for the group
		const new_group_element = document.createElement('div');

		// Set the class and ID of the new div element
		new_group_element.className = 'node-group';
		new_group_element.id = `node-group-${index + 1}`;

		// Set the inner HTML content of the new div element
		new_group_element.innerHTML = current_group_label;

		// Append the new div element to the container
		node_groups_container.appendChild(new_group_element);

		// Update the y-position and content of the new div element
		new_group_element.style.y = current_y_position;
		current_y_position += 206;
		new_group_element.innerHTML = `<span style="color: #ffc42c;">${current_group_label}</span><br><span style="color: #ffffff; margin-top: 9px; display: inline-block;">${label_data[3]} nodes</span>`;

		// Extract group properties from the label data and loop through each one
		const group_properties = label_data[0].split(',');
		for (const gp of group_properties) {
			const gp_id = `${label_data} -> ${gp}`;
			group_labels[gp] = gp;
			grouped_node_positions[gp_id] = label_data[2];
			node_radius_map[gp_id] = 12.5;
			node_colour_map[gp_id] = '#32CD32';
		}

		// Sanitize the group label by replacing certain substrings
		let sanitized_label = current_group_label.replace(': Supermodel', '');
		sanitized_label = sanitized_label.replace('1 & 2', '1&2');
		sanitized_label = sanitized_label.replace(/ • /g, '');
		sanitized_label = sanitized_label.replace(' entry', '');

		// Add the sanitized label to the nodes_data_array
		nodes_data_array.push({ id: current_group_label, label: sanitized_label, num_of_identicals: 150 });
		
		// Update the grouped_node_positions, node_radius_map, and node_colour_map with the new group label
		grouped_node_positions[current_group_label] = label_data[2];
		node_radius_map[current_group_label] = 20;
		node_colour_map[current_group_label] = 'none';

		// Increment the index (this is likely an unintended side effect and not necessary due to the forEach loop)
		index++;
	});
}

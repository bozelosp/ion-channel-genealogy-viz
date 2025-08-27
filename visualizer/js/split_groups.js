const MAX_GROUP_COUNT = 35;
const STARTING_Y_POSITION = 286;
const Y_POSITION_STEP = 206;
const LARGE_NODE_RADIUS = 20;
const SMALL_NODE_RADIUS = 12.5;

// Generates all unique combinations of elements from the given array.
// The combinations can be of any length from 1 to N, where N is the length of the array.
function get_combinations(arr) {
	let result = [];
	let f = function(prefix, arr) {
		for (let i = 0; i < arr.length; i++) {
			result.push([...prefix, arr[i]]);
			f([...prefix, arr[i]], arr.slice(i + 1));
		}
	};
	f([], arr);
	return result;
}

// Function to generate unique combinations of filters
function get_unique_filter_combinations(filter_state) {
	// Get active filters for Supermodel 1, Supermodel 2, and ICG entry
	const active_filters = ['Supermodel 1', 'Supermodel 2', 'ICG entry'].filter(f => filter_state[f]?.filter_value);

	// Get active ion classes based on the filter state
	const active_ion_classes = filter_state['All']?.filter_value ? ['All'] : ION_CLASS_NAMES.filter(c => filter_state[c]?.filter_value);

	// Initialize list to hold final unique combinations
	let unique_filter_combinations = [];

	// Generate all unique combinations of active filters
	const filter_combinations = get_combinations(active_filters);

	// Combine each filter combination with each ion class
	for (const ion_class of active_ion_classes) {
		for (const filter_combo of filter_combinations) {
			unique_filter_combinations.push([ion_class, ...filter_combo]);
		}
	}

	// Add single ion classes as their own group
	unique_filter_combinations = [...unique_filter_combinations, ...active_ion_classes.map(c => [c])];

	return unique_filter_combinations;
}

// Function to assign nodes to groups based on their properties
function assign_nodes_to_groups(nodes_data, unique_filter_combinations, filter_state) {

	// Initialize an object to store nodes grouped by filters
	const nodes_grouped_by_filter = {};

	// Loop through each unique filter combination
	unique_filter_combinations.forEach(filter_combination => {

		const group_key = filter_combination.join(',');
	
		// Loop through each node in the dataset
		nodes_data.forEach(node => {
			// Extract relevant properties from the node's 'original_model'

			let { ion_class, 'Supermodel 1': sm1, 'Supermodel 2': sm2, ICG: icg_entry } = node.original_model;

			let ion_class_match = filter_state['All']?.filter_value || filter_combination.includes(ion_class);

			let status = true ;

			if (ion_class_match) {

				if (filter_state['Supermodel 1'].filter_value === true) {

					if (filter_combination.includes('Supermodel 1')) {
						if (sm1!==true) {
							status = false ;
						}
					} else {
						if (sm1===true) {
							status = false ;
						}
					}
				}

				if (filter_state['Supermodel 2'].filter_value === true) {

					if (filter_combination.includes('Supermodel 2')) {
						if (sm2!==true) {
							status = false ;
						}
					} else {
						if (sm2===true) {
							status = false ;
						}
					}
				}

				if (filter_state['ICG entry'].filter_value === true) {

					if (filter_combination.includes('ICG entry')) {
						if (icg_entry!==true) {
							status = false ;
						}
					} else {
						if (icg_entry===true) {
							status = false ;
						}
					}
				}

				if (status===true) {

					nodes_grouped_by_filter[group_key] = nodes_grouped_by_filter[group_key] || []; // this makes sure that the key exists, if not it creates it and assigns an empty array
					nodes_grouped_by_filter[group_key].push(node); // this pushes the node id to the array of the corresponding key

				}
			
			}

		});
	
	});

	return nodes_grouped_by_filter;
}

function sort_and_position_groups(nodes_grouped_by_filter, fixed_location_circles) {

	console.log('nodes_grouped_by_filter',nodes_grouped_by_filter);

	let groups = [];
	for (const key in nodes_grouped_by_filter) {
		if (nodes_grouped_by_filter[key].length > 0) {
			const group_radius_values = nodes_grouped_by_filter[key].map(i => calculate_node_radius(i) );
			const total_group_radius = group_radius_values.reduce((a, b) => a + b, 0); // Added an initial value for the reduce function
			groups.push([key, Math.pow(total_group_radius, 0.8), nodes_grouped_by_filter[key]]);
		}
	}

	groups = groups.sort((a, b) => a[1] - b[1]); // Sort the groups by size in ascending order

	// Position the groups
	const group_circle_positions = position_circles(groups, fixed_location_circles); // Assuming positionCircles is another function you have

	const node_id_to_location = {};

	groups.forEach((group, index) => {

		const group_key = group[0];
		const group_radius = group[1];
		const group_nodes = group[2];

		const group_location = group_circle_positions[index];

		group_nodes.forEach(node => {
			node_id_to_location[node.id] = group_location;
		});

	});

	return node_id_to_location 

}

function position_circles(groups, fixed_location_circles) {

	const my_len = groups.length;
	const my_sum = groups.reduce((acc, curr) => acc + curr[1], 0);
	const normalized_d = groups.map(i => i[1] / my_sum);
	
	let my_min = Number.POSITIVE_INFINITY;
	let my_circle_positions;
	
	fixed_location_circles[my_len].forEach(i => {
		const my_score = i.reduce((acc, curr, index) => acc + Math.abs(curr[0] - normalized_d[index]), 0);
		if (my_score < my_min) {
		my_min = my_score;
		my_circle_positions = i;
		}
	});

	my_circle_positions = my_circle_positions.map( d => [avoid_borders(d[1]), avoid_borders(d[2])] );
		
	return my_circle_positions;
}

function avoid_borders(d){
	if (d>0.5){
		df=d-0.5;
		d-=-12*(Math.log(Math.sqrt(df))/Math.log(50))*Math.pow(df,2.5);
		return d;
	} else {
		df=0.5-d;
		d+=-12*(Math.log(Math.sqrt(df))/Math.log(50))*Math.pow(df,2.5);
		return d;
	}
}


function calculate_node_radius(d) {
	const my_radius = (0.2375 * Math.log(d.num_of_identicals) / Math.log(1.09) + 1.325 + 0.3925 * Math.log(d.num_of_identicals) / Math.log(1.35) + 3) /2 ;	
	return my_radius;
}


// Main function to group nodes based on selected filters
function group_nodes_by_selected_filters(nodes_data, filter_state, fixed_location_circles) {
	// Get unique filter combinations based on the current filter state
	const unique_filter_combinations = get_unique_filter_combinations(filter_state);

	// Assign nodes to groups based on these unique filter combinations
	const nodes_grouped_by_filter = assign_nodes_to_groups(nodes_data, unique_filter_combinations, filter_state);

	const number_of_groups = Object.keys(nodes_grouped_by_filter).length;

	if (number_of_groups > 1) {
		
		// Sort and position the groups (assuming you have a function for this)
		const node_id_to_location = sort_and_position_groups(nodes_grouped_by_filter, fixed_location_circles);

		return [nodes_grouped_by_filter, node_id_to_location];
	} else {

		return [nodes_grouped_by_filter, {}];
	
	}

}
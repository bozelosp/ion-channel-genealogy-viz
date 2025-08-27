// Constants for ion class names
const ION_CLASS_NAMES = ["K", "Na", "Ca", "Ih", "KCa", "Other"];

// Async function to initialize dashboard filters
function initialize_dashboard_filters(filter_state) {

	// Filter out checkboxes from filter state and extract relevant details
	const checkbox_data = Object.values(filter_state).filter(filter => filter.filter_type === "checkbox");
	const checkbox_ids = checkbox_data.map(filter => filter.filter_id);
	const checkbox_names = checkbox_data.map(filter => filter.filter_name);

	// Define IDs for numeric filter buttons
	const numeric_filter_button_ids = [
		"copies-number-decrease",
		"copies-number-increase",
		"similarity-score-decrease",
		"similarity-score-increase"
	];
  
	// Combine checkbox and numeric filter IDs
	const filter_ids = [...checkbox_ids, ...numeric_filter_button_ids];

	// Map numeric button IDs to filter names
	const numeric_button_id_to_name_map = {
		"copies-number-decrease": "# copies",
		"copies-number-increase": "# copies",
		"similarity-score-decrease": "similarity score",
		"similarity-score-increase": "similarity score"
	};

	// Map checkbox IDs to filter names
	const checkbox_id_to_name_map = Object.fromEntries(checkbox_data.map(filter => [filter.filter_id, filter.filter_name]));

	update_filter_display(checkbox_names,filter_state);

	return { 
		checkbox_names, 
		filter_state, 
		numeric_button_id_to_name_map, 
		checkbox_id_to_name_map, 
		filter_ids, 
		numeric_filter_button_ids 
	};
	
}

// Function to handle click events
async function handle_filter_click(event, item) {

	const { checkbox_names, filter_state, numeric_button_id_to_name_map, checkbox_id_to_name_map, filter_ids, numeric_filter_button_ids } = item;

	const clicked_id = event.target.id;

	const is_numeric_button = numeric_filter_button_ids.includes(clicked_id);

	// console.log(is_numeric_button);
	// console.log(clicked_id);
	// console.log(filter_state);

	const this_filter_name = is_numeric_button
		? numeric_button_id_to_name_map[clicked_id]
		: checkbox_id_to_name_map[clicked_id] ;

	// Branch logic based on the type of button clicked
	if (is_numeric_button) {
		handle_numeric_filter(clicked_id, this_filter_name);
	} else {
		handle_checkbox_filter(this_filter_name);
	}

	// // for each item in filter_state, log the key and the filter value
	// for (const [key, value] of Object.entries(filter_state)) {
	// 	console.log(`${key}: ${value.filter_value}`);
	// }

	// Update the DOM to reflect new filter state
	update_filter_display(checkbox_names,filter_state);

	// Function to handle checkbox filters
	function handle_checkbox_filter(this_filter_name) {

		// Toggle the current filter value
		filter_state[this_filter_name].filter_value =! filter_state[this_filter_name].filter_value;
		
		// Special case: if "All" is clicked, set all ion names to false
		if (this_filter_name === "All") {
			ION_CLASS_NAMES.forEach(function(d) {
				filter_state[d].filter_value = false;
			});
		} else if (ION_CLASS_NAMES.includes(this_filter_name)) {
			// If an individual ion name is clicked, set "All" to false
			filter_state["All"].filter_value = false;
		}

	}

	// Function to handle numeric filters
	function handle_numeric_filter(clicked_id, this_filter_name) {
	
		// Define handlers for each numeric filter button
		const handlers = {
			'copies-number-decrease': () => Math.max(filter_state['# copies'].min, filter_state['# copies'].filter_value - 1),
			'copies-number-increase': () => Math.min(filter_state['# copies'].max, filter_state['# copies'].filter_value + 1),
			'similarity-score-decrease': () => Math.max(75, filter_state['similarity score'].filter_value - 1),
			'similarity-score-increase': () => Math.min(100, filter_state['similarity score'].filter_value + 1)
		};

		// Execute the appropriate handler based on the clicked button
		if (clicked_id in handlers) {
			filter_state[this_filter_name].filter_value = handlers[clicked_id]();
		}
	}

}

// Function to update the DOM based on the new filter state
function update_filter_display(checkbox_names,filter_state) {
	
	// Update checkbox elements
	checkbox_names.forEach(function(d) {
		let elem = document.getElementById(filter_state[d].filter_id);
		elem.classList.toggle("active", filter_state[d].filter_value);
	});

	// Update numeric fields
	document.getElementById('copies-number').value = filter_state['# copies'].filter_value;
	document.getElementById('similarity-score-number').value = filter_state['similarity score'].filter_value;
}
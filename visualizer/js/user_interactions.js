async function setup_click_events_for_nodes(links_data, node, network_data, fetched_files) {

	d3.selectAll("circle").on("click", function() {
		const current_clicked_node = d3.select(this);

		[ source_node_ids, target_node_ids ] = [ [], [] ];

		// rmeove selected-node class from all nodes
		d3.selectAll(".selected-node").classed("selected-node", false);
		// remove source-node class from all nodes
		d3.selectAll(".source-node").classed("source-node", false);

		console.log("current_clicked_node: ", current_clicked_node);
		console.log("event: ", d3.event);

		if (d3.event.shiftKey) {

			let subgraph_node_ids = [];
			current_clicked_node.each(node_data => subgraph_node_ids.push(node_data.id));

			let copied_links_data = Array.from(links_data);

			let previous_subgraph_size = -1;
			let current_subgraph_size = 1;

			while (previous_subgraph_size !== current_subgraph_size) {
			
				previous_subgraph_size = current_subgraph_size;
				extend_subgraph_by_connected_nodes(copied_links_data, subgraph_node_ids);
				subgraph_node_ids = [...new Set(subgraph_node_ids)];  // Remove duplicates here
				current_subgraph_size = subgraph_node_ids.length;
			}

			//let modeldb_time_prefixes = subgraph_node_ids.map(id => parseInt(/^(\d+)/g.exec(id)[1]));
			// we can actually extract the time info from hear:
			// {
			// 	"original_model":
			// 	{
			// 		"mod_filepath": "modelDB_unzipped_models/258946/Alpha5_NMDA_single_comp/nmdaSyn.mod",
			// 		"mod_filename": "nmdaSyn.mod",
			// 		"unique_modelDB_mod_id": "258946_nmdaSyn-ID-1",
			// 		"modelDB_dir": "258946",
			// 		"ICG": false,
			// 		"ICG_entries": null,
			// 		"ion_class": null,
			// 		"Supermodel 1": true,
			// 		"Supermodel 2": true,
			// 		"Year": 2018
			// 	}, 
			
			// Get all nodes in the subgraph
			const nodes_in_subgraph = node.filter(node_data => subgraph_node_ids.includes(node_data.id));
			nodes_in_subgraph.classed("selected-node", true);

			// Extract the year information for each node in the subgraph
			const years_in_subgraph = nodes_in_subgraph.data().map(node_data => node_data.original_model.Year ? node_data.original_model.Year : 2023);

			// Find the minimum (earliest) year
			const earliest_year = Math.min(...years_in_subgraph);
			
			// Filter node by the earliest year
			const original_node_in_subgraph = nodes_in_subgraph.filter(node_data => node_data.original_model.Year === earliest_year);

			// if more than one node has the same earliest year, then choose the one with the smallest modelDB_dir parsed as an integer
			if (original_node_in_subgraph.size() > 1) {
				const modelDB_dirs = original_node_in_subgraph.data().map(node_data => parseInt(node_data.original_model.modelDB_dir));
				const smallest_modelDB_dir = Math.min(...modelDB_dirs);
				original_node_in_subgraph.filter(node_data => parseInt(node_data.original_model.modelDB_dir) === smallest_modelDB_dir);
			}

			// Class the earliest node as the source node
			original_node_in_subgraph.classed("source-node", true);

			// return the original node id(s) and the rest of the subgraph node ids
			source_node_ids = original_node_in_subgraph.data().map(node_data => node_data.id);
			target_node_ids = subgraph_node_ids.filter(node_id => !source_node_ids.includes(node_id));

			// console.log("source_node_ids: ", source_node_ids);
			// console.log("target_node_ids: ", target_node_ids);

			const fetch_these_node_ids = [...source_node_ids, ...target_node_ids];

			fetch_source_code(fetch_these_node_ids, network_data, fetched_files) ;
			
			setup_keyboard_interactions(source_node_ids, target_node_ids, fetched_files);

		} else if (d3.event.ctrlKey || d3.event.metaKey) {
			current_clicked_node.classed("source-node", !current_clicked_node.classed("source-node"));
		} else {
			current_clicked_node.classed("selected-node", !current_clicked_node.classed("selected-node"));
		}
	});

}
function setup_keyboard_interactions(source_node_ids, target_node_ids, fetched_files) {
    document.addEventListener("keyup", function(event) {
        if (event.code === 'KeyD') {
            handle_diff(source_node_ids, target_node_ids, fetched_files);

            // Adding an inner event listener for arrow keys
            document.addEventListener('keyup', function(event) {
                switch (event.code) {
                    case 'ArrowLeft':
                        navigate_diff_boxes('left');
                        break;
                    case 'ArrowRight':
                        navigate_diff_boxes('right');
                        break;
                    default:
                        break;
                }
            });
        }
    });
}


function extend_subgraph_by_connected_nodes(link_data_array, subgraph_node_ids) {
	link_data_array.forEach(link => {
		if (subgraph_node_ids.includes(link.source.id) || subgraph_node_ids.includes(link.target.id)) {
			subgraph_node_ids.push(link.source.id);
			subgraph_node_ids.push(link.target.id);
		}
	});
	subgraph_node_ids = [...new Set(subgraph_node_ids)];  // Remove duplicates
}

function apply_css_classes_to_nodes() {
	// Moved the styles to CSS, so no longer need to set them here
	d3.selectAll(".source-node, .selected-node")
		.transition()
		.duration(212.5);
}






function setup_drag_handler(node) {
	let drag_handler = d3.drag()
		.on("start", drag_started)
		.on("drag", dragging)
		.on("end", drag_ended);
	drag_handler(node);
}

function setup_mouse_events() {
	d3.selectAll("circle")
		.on("mouseover", function(d){

			let str = "<span style=\"color: #ffc42c; font-size: 1.18em;\">Details</span><br>" ;

				str += "<span style=\"color: #ffc42c; line-height: 40px;\">ICG id: <\/span>" + d.id + "<br>"
					+ "<span style=\"color: #ffc42c; line-height: 28px;\">Self + identicals: <\/span>" + d.num_of_identicals + "<br>";

			d3.select("#mouseover-details-box")
				.html(str);

			str="";
			d.list_of_identicals.forEach(function(x){
					str = str + "<span style=\"line-height: 26px\">" + x + "<span><br>" ;
			});

			d3.select('#siblings-family-box')
				.html(str);

		})
		.on("mouseout", function(d){});

		recolour_nodes();
	
}

// function setup_click_events() {

//     d3.selectAll("circle").on("click", function() {

// 		if (d3.event.shiftKey) {
		
// 			let clicked_node_subgraph = [];
// 			d3.select(this).each(d => clicked_node_subgraph.push(d.id));
	
// 			let my_links_data = Array.from(links_data); // Assuming links_data is a Map or Set
	
// 			let previous_arr_length = -1;
// 			let current_arr_length = 1;
	
// 			while (previous_arr_length !== current_arr_length) {
// 				previous_arr_length = current_arr_length;
// 				retrieve_subgraph_node_ids(my_links_data);
// 				current_arr_length = clicked_node_subgraph.length;
// 			}
	
// 			let modeldb_id_time_arr = clicked_node_subgraph.map(d => parseInt(/^(\d+)/g.exec(d)[1]));
// 			let original_node_modeldb_id = Math.min(...modeldb_id_time_arr);
	
// 			let selected_subgraph = node.filter(d => clicked_node_subgraph.includes(d.id));
// 			selected_subgraph.attr("class", "selected-node");
	
// 			let original_node = selected_subgraph.filter(d => d.id.includes(`${original_node_modeldb_id}-`));
// 			original_node.attr("class", "source-node");
	
// 		}			

// 		else if (d3.event.ctrlKey || d3.event.metaKey) {
// 			if (this.getAttribute("class") !== "source-node") {
// 				d3.select(this).attr("class", "source-node");
// 			}
// 		} else {
// 			let currentClass = this.getAttribute("class");
// 			if (currentClass === "source-node") {
// 				d3.select(this).attr("class", "deselected-node");
// 			} else if (currentClass === "selected-node") {
// 				d3.select(this).attr("class", "deselected-node");
// 			} else {
// 				d3.select(this).attr("class", "selected-node");
// 			}
// 		}
// 	});
// }

// function recolour_nodes() {

// 	d3.selectAll(".source-node")
// 		.transition()
// 		.duration(212.5)
// 		.style("fill", "#ffd700")
// 		.style("stroke", "#215885")
// 		.style("stroke-width","2")
// 		.attr("class","source-node");

// 	d3.selectAll(".selected-node")
// 		.transition()
// 		.duration(212.5)
// 		.style("fill", "#00BFFF")
// 		.style("stroke", "#215885")
// 		.style("stroke-width","2");
	
// 	d3.selectAll(".deselected-node")
// 		.transition()
// 		.duration(212.5)
// 		.style("fill", "#00BFFF")
// 		.style("stroke", "#aaa")
// 		.style("stroke-width","1");
// }

// function retrieve_subgraph_node_ids(my_arr,clicked_node_subgraph) {
// 	my_arr.forEach(d => {
// 		if (clicked_node_subgraph.includes(d.source.id) || clicked_node_subgraph.includes(d.target.id)) {
// 			clicked_node_subgraph.push(d.source.id);
// 			clicked_node_subgraph.push(d.target.id);
// 		}
// 	});
// 	clicked_node_subgraph = clicked_node_subgraph.filter(onlyUnique);
// }
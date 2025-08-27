// A function to handle the "Enter" key action
function handleEnterKey() {
    let my_query = document.getElementById('my-input').value.split(' ');

		// Here we assume that 'node' is a global variable or in scope
		let queried_nodes = node.filter(function(d) {
			let my_status = false;
			let mod_files_represented = d.list_of_identicals.concat([d.id]).filter(onlyUnique);

			mod_files_represented.forEach(function(x) {
				let ttt = x.match(/^(\d+)/)[1];
				my_query.forEach(function(r) {
					if (r === ttt) {
						my_status = true;
					}
				});
			});
			return my_status;
		});

		queried_nodes.forEach(node => {
			node.className = 'selected-node';
		});

		let clicked_node_subgraph = queried_nodes.map(node => node.id);

		// Assuming my_links_data is already defined in the code
		let my_links_data = Object.values(links_data);

		let previous_arr_length = -1;
		let current_arr_length = 1;

		while (previous_arr_length !== current_arr_length) {
			previous_arr_length = current_arr_length;
			// Assuming retrieve_subgraph_node_ids() is a function defined elsewhere in your code
			retrieve_subgraph_node_ids(my_links_data);
			current_arr_length = clicked_node_subgraph.length;
		}

		let modeldb_id_time_arr = clicked_node_subgraph.map(d => parseInt(d.match(/^(\d+)/)[1]));

		let original_node_modeldb_id = Math.min(...modeldb_id_time_arr);

		// Assuming selected_subgraph is an array or similar structure
		let selected_subgraph = node.filter(d => clicked_node_subgraph.includes(d.id));

		selected_subgraph.forEach(node => {
			node.className = 'selected-node';
		});

		let original_node = selected_subgraph.filter(d => d.id.includes(`${original_node_modeldb_id}-`));

		original_node.forEach(node => {
			node.className = 'source-node';
		});

		// Assuming recolour_nodes() is a function defined elsewhere in your code
		recolour_nodes();
}

// A function to handle the "Backspace" or "Escape" key action
function handleBackspaceOrEscape() {
    d3.selectAll(".selected-node")
			.transition()
			.duration(212.5)
			.style("fill", "#00BFFF")
			.style("stroke", "#aaa")
			.style("stroke-width","1")
			.attr("class","deselected-node");

		d3.selectAll(".source-node")
			.transition()
			.duration(212.5)
			.style("fill", "#00BFFF")
			.style("stroke", "#aaa")
			.style("stroke-width","1")
			.attr("class","deselected-node");

		difflib_container_box = document.getElementById('difflib-container') ;

		for(jj=1; jj < 30; jj++) {
			select_this_difflib_box = 'difflib-box-' + jj.toString() ;
			d3.select("#" + select_this_difflib_box).remove() ;	
		}

		jj_max = null ;

		document.getElementById('diff-source-label').innerHTML = '';
		document.getElementById('diff-target-label').innerHTML = '';
		document.getElementById('diff-id-label').innerHTML = '';
}

// A function to handle the "D" or "d" key action
function handle_D_Key() {
    try {
		const target_nodes = Array.from(document.querySelectorAll('.selected-node')).map(d => d.id);
		const source_nodes = Array.from(document.querySelectorAll('.source-node')).map(d => d.id);

		const difflib_container_box = document.getElementById('difflib-container');

		for (let jj = 1; jj < 30; jj++) {
			const select_this_difflib_box = 'difflib-box-' + jj.toString();
			const element_to_remove = document.getElementById(select_this_difflib_box);
			if (element_to_remove) element_to_remove.remove();
		}

		let jj = 1;
		let difflib_table_ids = [];
		const my_reg_exp = /(difflib_chg.*?top)/g;

		source_nodes.forEach(mySourceNode => {
			target_nodes.forEach(myTargetNode => {
				const select_this_difflib_box = 'difflib-box-' + jj.toString();

				let new_element = document.createElement('div');
				new_element.className = `difflib-box ${mySourceNode} ${myTargetNode}`;
				new_element.id = select_this_difflib_box;

				const difflibF1F2Path = 'https://.../' + mySourceNode + '/' + mySourceNode + '_' + myTargetNode + '.html';

				fetch(difflibF1F2Path)
					.then(response => response.text())
					.then(myHtml => {
						const match = my_reg_exp.exec(myHtml);
						const my_match = match[1];
						difflib_table_ids.push(my_match);

						new_element.innerHTML = myHtml;
						const table = new_element.querySelector("#" + my_match);
						new_element.innerHTML = table.outerHTML;

						const elems = new_element.querySelectorAll('td:not(.diff_header)');
						elems.forEach(elem => {
							elem.style.width = diffColumnWidth;
						});
					});

				difflib_container_box.appendChild(new_element);

				if (jj > 1) {
					document.querySelector('#' + select_this_difflib_box).style.visibility = 'hidden';
				}

				jj++;
			});
		});

		const jj_max = jj - 1;
		const select_this_difflib_box = 'difflib-box-1';
		const cls_nm = document.getElementById(select_this_difflib_box).className;
		const cls_nm_arr = cls_nm.split(" ");

		document.getElementById('diff-source-label').text_content = cls_nm_arr[1];
		document.getElementById('diff-target-label').text_content = cls_nm_arr[2];
		document.getElementById('diff-id-label').text_content = `1 / ${jj_max}`;
	} catch (error) {
		console.log('check the function that gets the node ids of the selected source and target .mod files to show difflib for');
	}
}

function handle_arrow_navigation(direction, jj_max, diff_where_ami) {
    let boundary_condition = false;
    let step = 0;
    
    if (direction === "left") {
        boundary_condition = (jj_max !== null && diff_where_ami !== 1);
        step = -1;
    } else if (direction === "right") {
        boundary_condition = (jj_max !== null && diff_where_ami < jj_max);
        step = 1;
    }
    
    if (boundary_condition) {
        const hide_element_id = `difflib-box-${diff_where_ami}`;
        document.getElementById(hide_element_id).style.visibility = 'hidden';

        diff_where_ami += step;

        const show_element_id = `difflib-box-${diff_where_ami}`;
        document.getElementById(show_element_id).style.visibility = 'visible';

        const cls_nm = document.getElementById(show_element_id).className;
        const cls_nm_arr = cls_nm.split(" ");

        document.getElementById('diff-source-label').text_content = cls_nm_arr[1];
        document.getElementById('diff-target-label').text_content = cls_nm_arr[2];
        document.getElementById('diff-id-label').text_content = `${diff_where_ami} / ${jj_max}`;
    }
}

// Map keys to their corresponding actions
const keyActions = {
    'Enter': handleEnterKey,
    'Backspace': handleBackspaceOrEscape,
    'Escape': handleBackspaceOrEscape,
    'D': handle_D_Key,
    'd': handle_D_Key
};

// In your event listener:
document.addEventListener("keyup", function(event) {
    const action = keyActions[event.key];
    if (action) {
        action();
    } else if (event.key === "ArrowLeft") {
        handle_arrow_navigation("left", jj_max, diff_where_ami);
    } else if (event.key === "ArrowRight") {
        handle_arrow_navigation("right", jj_max, diff_where_ami);
    }
});
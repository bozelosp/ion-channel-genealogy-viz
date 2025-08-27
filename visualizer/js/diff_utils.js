let max_diff_index = 0; // Global counter for the max number of diff boxes
let current_diff_index = 1; // Global counter for the current diff box being displayed

// Function to handle diff creation and display
async function handle_diff(source_node_ids, target_node_ids, fetched_files) {

	const diff_column_width = Math.round((document.documentElement.clientWidth - 755) / 2);
	const diff_box_width = document.documentElement.clientWidth - 680;

	console.log('diff_column_width', diff_column_width);

	const diff_container = document.getElementById('difflib-container');
	diff_container.innerHTML = ''; // Clear any existing diff boxes
	max_diff_index = 0; // Reset max_diff_index

	for (const source_id of source_node_ids) {
		const source_code = fetched_files[source_id]?.source_code;
		if (!source_code) {
			console.warn(`Source code for node ${source_id} is not fetched.`);
			continue;
		}

		for (const target_id of target_node_ids) {
			const target_code = fetched_files[target_id]?.source_code;
			if (!target_code) {
				console.warn(`Target code for node ${target_id} is not fetched.`);
				continue;
			}

			const diff_table_string = await request_diff_html(source_code, target_code);

			const diff_box_element = document.createElement('div');
			// add class to the diff-box
			diff_box_element.classList.add('difflib-box');
			diff_box_element.style.width = `${diff_box_width}px`;
			diff_box_element.innerHTML = diff_table_string;

			const myRegexp = /(difflib_chg.*?top)/g;
			const match = myRegexp.exec(diff_table_string) ;
			const my_match = match[1] ;
			const table=diff_box_element.querySelector("#"+my_match) ;
			diff_box_element.innerHTML = table.outerHTML ;
			d3.selectAll('.diff_next').remove();

			// Select all <td> elements with the 'nowrap' attribute
			const tdElements = document.querySelectorAll('td:not(.diff_header)');
			tdElements.forEach(elem => elem.removeAttribute('nowrap'));
			tdElements.forEach(elem => elem.style.width = `${diff_column_width}px`);

			diff_box_element.style.visibility = 'hidden'; // Hide by default
			diff_box_element.id = `difflib-box-${++max_diff_index}`; // Assign an ID
			diff_container.appendChild(diff_box_element);
		}
	}

	// Show the first diff box
	document.getElementById('difflib-box-1').style.visibility = 'visible';
	current_diff_index = 1; // Reset to the first index

	// Update the diff labels
	update_diff_labels();
}

// Function to request diff HTML from AWS Lambda
async function request_diff_html(source_code, target_code) {
	try {
		const diff_url = `https://el3l8it47l.execute-api.eu-west-1.amazonaws.com/production/generateHtmlDiff?string1=${encodeURIComponent(source_code)}&string2=${encodeURIComponent(target_code)}`;
		const response = await fetch(diff_url);
		if (!response.ok) {
			throw new Error('Network response was not ok');
		}
		const prepend_html_text = prepend_html();
		const diff_html_text = await response.text();

		return prepend_html_text + diff_html_text;
	} catch (error) {
		console.error('Fetch error:', error);
	}
}

// Update labels showing current diff box
function update_diff_labels() {
	document.getElementById('diff-id-label').innerHTML = `${current_diff_index} / ${max_diff_index}`;
}

// Function to navigate between diff boxes
function navigate_diff_boxes(direction) {
	// Hide the current diff box
	document.getElementById(`difflib-box-${current_diff_index}`).style.visibility = 'hidden';
	
	// Update the current index based on the direction
	if (direction === 'left' && current_diff_index > 1) {
		current_diff_index--;
	} else if (direction === 'right' && current_diff_index < max_diff_index) {
		current_diff_index++;
	}

	// Show the new current diff box
	document.getElementById(`difflib-box-${current_diff_index}`).style.visibility = 'visible';
	
	// Update the diff labels
	update_diff_labels();
}

function prepend_html() {

	const html_text = `<!DOCTYPE html PUBLIC "-//W3C//DTD XHTML 1.0 Transitional//EN"
	"http://www.w3.org/TR/xhtml1/DTD/xhtml1-transitional.dtd">
	
	<html>
	
	<head>
			<meta http-equiv="Content-Type"
						content="text/html; charset=utf-8" />
			<title></title>
			<style type="text/css">
					table.diff {font-family:Courier; border:medium;}
					.diff_header {background-color:#e0e0e0}
					td.diff_header {text-align:right}
					.diff_next {background-color:#c0c0c0}
					.diff_add {background-color:#aaffaa}
					.diff_chg {background-color:#ffff77}
					.diff_sub {background-color:#ffaaaa}
			</style>
	</head>
	
	<body>'''
	
		my_tail='''    <table class="diff" summary="Legends">
					<tr> <th colspan="2"> Legends </th> </tr>
					<tr> <td> <table border="" summary="Colors">
												<tr><th> Colors </th> </tr>
												<tr><td class="diff_add">&nbsp;Added&nbsp;</td></tr>
												<tr><td class="diff_chg">Changed</td> </tr>
												<tr><td class="diff_sub">Deleted</td> </tr>
										</table></td>
							 <td> <table border="" summary="Links">
												<tr><th colspan="2"> Links </th> </tr>
												<tr><td>(f)irst change</td> </tr>
												<tr><td>(n)ext change</td> </tr>
												<tr><td>(t)op</td> </tr>
										</table></td> </tr>
			</table>
	</body>
	
	</html>` ;

	return html_text;
}
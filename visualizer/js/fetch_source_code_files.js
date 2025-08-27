// Function to handle the mouseover event
async function fetch_source_code(node_ids_array, network_data, fetched_files) {

    const base_url = "http://ion-channels.s3-website-eu-west-1.amazonaws.com/static/modelDB/";

    for (const node_id of node_ids_array) {
        // Find the node data corresponding to the current node ID
        const node_data = network_data.nodes.find(node => node.id === node_id);
        if (!node_data) {
            console.warn(`Node with ID ${node_id} not found.`);
            continue;
        }

        const unique_id = node_data.original_model.unique_modelDB_mod_id;

        // Check if the files for this node have already been fetched
        if (!fetched_files[node_id]) {

            // File paths based on your AWS S3 bucket structure
            const source_code_path = `${base_url}source_code/${unique_id}.mod`;
            const stripped_comments_path = `${base_url}stripped_comments/${unique_id}.mod`;
            const stripped_whitespace_path = `${base_url}stripped_comments_whitespaces/${unique_id}.mod`;

            // Fetch the files and save them in the object
            await Promise.all([
                fetch(source_code_path).then(response => response.text()),
                fetch(stripped_comments_path).then(response => response.text()),
                fetch(stripped_whitespace_path).then(response => response.text())
            ]).then(([source_code, stripped_comments, stripped_whitespace]) => {
                fetched_files[node_id] = {
                    source_code: source_code
                    // source_code_stripped_comments: stripped_comments,
                    // source_code_stripped_whitespace: stripped_whitespace
                };
            }).catch(error => {
                console.error(`Error fetching files for node ${node_id}:`, error);
            });
        }
    }

    return fetched_files;

}
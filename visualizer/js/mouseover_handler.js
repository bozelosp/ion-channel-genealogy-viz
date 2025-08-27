// mouseover_handler.js

function handle_node_mouseover(node, nodes_data) {
    node.on("mouseover", function(event, d) {
        const modelDB_id = nodes_data[d].original_model.unique_modelDB_mod_id;
        const num_of_identicals = nodes_data[d].num_of_identicals;

        const detailsStr = `
            <span style="color: #ffc42c; font-size: 1.18em;">Details</span><br>
            <span style="color: #ffc42c; line-height: 40px;">ModelDB ID: </span>${modelDB_id}<br>
            <span style="color: #ffc42c; line-height: 28px;">Self + identicals: </span>${num_of_identicals}<br>
        `;
    
        d3.select("#mouseover-details-box").html(detailsStr);

        const siblingsStr = nodes_data[d].identical_models.map(model => `
            <span style="line-height: 26px">${model.unique_modelDB_mod_id}</span><br>
        `).join('');

        d3.select('#siblings-family-box').html(siblingsStr);
    });
}
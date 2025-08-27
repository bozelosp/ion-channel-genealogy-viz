// Initialize variables related to graph dimensions and layout

const this_width = document.documentElement.clientWidth;
const this_height = document.documentElement.clientHeight -40;

const graph_width = this_width - 510;
const graph_height = this_height - 9;
// Set SVG dimensions for the graph
document.getElementById('my_graph').setAttribute('width', graph_width);
document.getElementById('my_graph').setAttribute('height', graph_height);

// Set dimensions for gradient rectangles
document.getElementById('gradient_rect_left').setAttribute('width', this_width * 0.1);
document.getElementById('gradient_rect_left').setAttribute('height', this_height);

// Set heights for different layout components
document.getElementById('gradient_rect_right').setAttribute('width', this_width * 0.1);
document.getElementById('gradient_rect_right').setAttribute('height', this_height);

const width = document.documentElement.clientWidth;
const height = document.documentElement.clientHeight;

const diff_box_center = 296 + (width - 680) / 2 ;
const diff_column_width = Math.round((width - 755) / 2);
const label_offset = (33.5 + diff_column_width) / 2;

// Set widths for the diff utility and search bar
document.getElementById('node-groups-box').style.height = `${height - 318}px`;
document.getElementById('siblings-family-box').style.height = `${height - 259}px`;
document.getElementById('gradient_rect_bottom_right').style.height = `${height - 467}px`;

// Set widths for the diff utility and search bar
document.getElementById('difflib-container').style.width = `${width - 680}px`;
document.getElementById('search-bar').style.width = `${width - 680}px`;
document.getElementById('my-input').style.width = `${width - 700}px`;

// Set positions for the labels in the diff utility
document.getElementById('diff-source-label').style.left = `${diff_box_center - label_offset - 33.5-10.5}px`;
document.getElementById('diff-target-label').style.left = `${diff_box_center + label_offset - 33.5-7}px`;
document.getElementById('diff-id-label').style.left = `${width - 680 + 296 - 37}px`;
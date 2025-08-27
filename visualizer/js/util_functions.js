// Separate JavaScript file, for example 'helperFunctions.js'

var radius = 3.25;

//d3.select(window).on('resize.updatesvg', updateWindow);

var hexDigits = new Array("0","1","2","3","4","5","6","7","8","9","A","B","C","D","E","F"); 



function dragStarted(simulation, d) {
	if (!d3.event.active) simulation.alphaTarget(0.3).restart();
	d.fx = d.x;
	d.fy = d.y;
}

function dragging(d) {
	d.fx = d3.event.x;
	d.fy = d3.event.y;
}

function dragEnded(simulation, d) {
	if (!d3.event.active) simulation.alphaTarget(0);
	d.fx = null;
	d.fy = null;
}

function zoomed(g) {
	g.attr("transform", d3.event.transform);
}

function getKeyByValue(object, value) {
	return "." + Object.keys(object).find(key => object[key] === value);
}

function onlyUnique(value, index, self) {
	return self.indexOf(value) === index;
}

function rgb2hex(rgb) {
	rgb = rgb.match(/^rgb\((\d+),\s*(\d+),\s*(\d+)\)$/);
	return "#" + hex(rgb[1]) + hex(rgb[2]) + hex(rgb[3]);
}

function hex(x) {
	return isNaN(x) ? "00" : hexDigits[(x - x % 16) / 16] + hexDigits[x % 16];
}
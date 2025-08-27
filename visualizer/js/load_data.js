async function loadData() {
	
	const urls = [
		'data/filter_state.json',
		'data/fixed_location_circles.json',
		'data/network_data.json'
	];
	
	const results = [];
  
	const fetch_promises = urls.map(async (url, index) => {
	  try {
		const response = await fetch(url);
		if (!response.ok) {
		  throw new Error(`Failed to fetch ${url}: ${response.statusText}`);
		}
		results[index] = await response.json();
	  } catch (error) {
		console.error(`Error loading JSON data from ${url}:`, error);
	  }
	});
  
	await Promise.all(fetch_promises);
	
	return results; // [filter_state, fixed_location_circles, network_data]
}
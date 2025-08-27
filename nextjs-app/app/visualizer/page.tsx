'use client';

import { useEffect, useRef, useState } from 'react';
import * as d3 from 'd3';

interface Node {
  id: string;
  name: string;
  ion_class?: string;
  supermodel?: number;
  icg?: boolean;
  num_of_identicals?: number;
  original_model?: {
    ion_class?: string;
    ICG?: boolean;
  };
  x?: number;
  y?: number;
  fx?: number | null;
  fy?: number | null;
}

interface Link {
  source: string | Node;
  target: string | Node;
  value?: number;
  weight?: number;
}

interface NetworkData {
  nodes: Node[];
  links: Link[];
}

// Calculate node radius based on num_of_identicals (from original implementation)
function calculateNodeRadius(node: Node): number {
  const numIdenticals = node.num_of_identicals || 1;
  const radius = (0.2375 * Math.log(numIdenticals) / Math.log(1.09) + 1.325 + 
                  0.3925 * Math.log(numIdenticals) / Math.log(1.35) + 3) / 2;
  return radius;
}

// Group splitting constants
const ION_CLASS_NAMES = ['K', 'Na', 'Ca', 'IH', 'KCa', 'Other'];
const FORCE_STRENGTH = 0.1675;

// Get all combinations of array elements
function getCombinations(arr: string[]): string[][] {
  const result: string[][] = [];
  const f = (prefix: string[], arr: string[]) => {
    for (let i = 0; i < arr.length; i++) {
      result.push([...prefix, arr[i]]);
      f([...prefix, arr[i]], arr.slice(i + 1));
    }
  };
  f([], arr);
  return result;
}

// Get unique filter combinations based on active filters
function getUniqueFilterCombinations(ionClass: string, showICG: boolean, supermodel1: boolean, supermodel2: boolean): string[][] {
  // Get active filters
  const activeFilters: string[] = [];
  if (supermodel1) activeFilters.push('Supermodel 1');
  if (supermodel2) activeFilters.push('Supermodel 2');
  if (showICG) activeFilters.push('ICG entry');

  // Get active ion classes
  const activeIonClasses = ionClass === 'all' ? ['All'] : [ionClass];

  // Initialize list for unique combinations
  let uniqueFilterCombinations: string[][] = [];

  // Generate all unique combinations of active filters
  const filterCombinations = getCombinations(activeFilters);

  // Combine each filter combination with each ion class
  for (const ionCls of activeIonClasses) {
    for (const filterCombo of filterCombinations) {
      uniqueFilterCombinations.push([ionCls, ...filterCombo]);
    }
  }

  // Add single ion classes as their own group
  uniqueFilterCombinations = [...uniqueFilterCombinations, ...activeIonClasses.map(c => [c])];

  return uniqueFilterCombinations;
}

// Assign nodes to groups based on their properties
function assignNodesToGroups(
  nodes: Node[], 
  uniqueFilterCombinations: string[][], 
  ionClass: string,
  showICG: boolean,
  supermodel1: boolean,
  supermodel2: boolean
): { [key: string]: Node[] } {
  const nodesGroupedByFilter: { [key: string]: Node[] } = {};

  uniqueFilterCombinations.forEach(filterCombination => {
    const groupKey = filterCombination.join(',');

    nodes.forEach(node => {
      const nodeIonClass = node.ion_class || node.original_model?.ion_class;
      const nodeSm1 = node.supermodel === 1;
      const nodeSm2 = node.supermodel === 2;
      const nodeIcg = node.icg || node.original_model?.ICG || false;

      const ionClassMatch = ionClass === 'all' || filterCombination.includes(nodeIonClass || '') || filterCombination.includes('All');

      let status = true;

      if (ionClassMatch) {
        // Check Supermodel 1
        if (supermodel1) {
          if (filterCombination.includes('Supermodel 1')) {
            if (!nodeSm1) status = false;
          } else {
            if (nodeSm1) status = false;
          }
        }

        // Check Supermodel 2
        if (supermodel2) {
          if (filterCombination.includes('Supermodel 2')) {
            if (!nodeSm2) status = false;
          } else {
            if (nodeSm2) status = false;
          }
        }

        // Check ICG entry
        if (showICG) {
          if (filterCombination.includes('ICG entry')) {
            if (!nodeIcg) status = false;
          } else {
            if (nodeIcg) status = false;
          }
        }

        if (status) {
          if (!nodesGroupedByFilter[groupKey]) {
            nodesGroupedByFilter[groupKey] = [];
          }
          nodesGroupedByFilter[groupKey].push(node);
        }
      }
    });
  });

  return nodesGroupedByFilter;
}

// Avoid borders function from original
function avoidBorders(d: number): number {
  if (d > 0.5) {
    const df = d - 0.5;
    d -= -12 * (Math.log(Math.sqrt(df)) / Math.log(50)) * Math.pow(df, 2.5);
    return d;
  } else {
    const df = 0.5 - d;
    d += -12 * (Math.log(Math.sqrt(df)) / Math.log(50)) * Math.pow(df, 2.5);
    return d;
  }
}

// Position circles based on groups
function positionCircles(groups: any[], fixedLocationCircles: any): [number, number][] {
  const myLen = groups.length;
  const mySum = groups.reduce((acc, curr) => acc + curr[1], 0);
  const normalizedD = groups.map(i => i[1] / mySum);
  
  let myMin = Number.POSITIVE_INFINITY;
  let myCirclePositions;
  
  if (fixedLocationCircles && fixedLocationCircles[myLen]) {
    fixedLocationCircles[myLen].forEach((i: any) => {
      const myScore = i.reduce((acc: number, curr: any, index: number) => 
        acc + Math.abs(curr[0] - normalizedD[index]), 0);
      if (myScore < myMin) {
        myMin = myScore;
        myCirclePositions = i;
      }
    });
  }

  if (!myCirclePositions) {
    // Fallback positioning if no fixed locations available
    return groups.map((_, index) => [0.5, 0.5]);
  }

  return myCirclePositions.map((d: any) => [avoidBorders(d[1]), avoidBorders(d[2])]);
}

// Sort and position groups
function sortAndPositionGroups(
  nodesGroupedByFilter: { [key: string]: Node[] },
  fixedLocationCircles: any
): { [key: string]: [number, number] } {
  const groups: any[] = [];
  
  for (const key in nodesGroupedByFilter) {
    if (nodesGroupedByFilter[key].length > 0) {
      const groupRadiusValues = nodesGroupedByFilter[key].map(node => calculateNodeRadius(node));
      const totalGroupRadius = groupRadiusValues.reduce((a, b) => a + b, 0);
      groups.push([key, Math.pow(totalGroupRadius, 0.8), nodesGroupedByFilter[key]]);
    }
  }

  // Sort groups by size
  groups.sort((a, b) => a[1] - b[1]);

  // Position the groups
  const groupCirclePositions = positionCircles(groups, fixedLocationCircles);

  const nodeIdToLocation: { [key: string]: [number, number] } = {};

  groups.forEach((group, index) => {
    const groupNodes = group[2];
    const groupLocation = groupCirclePositions[index];

    groupNodes.forEach((node: Node) => {
      nodeIdToLocation[node.id] = groupLocation;
    });
  });

  return nodeIdToLocation;
}

export default function Visualizer() {
  const svgRef = useRef<SVGSVGElement>(null);
  const [networkData, setNetworkData] = useState<NetworkData | null>(null);
  const [fixedLocationCircles, setFixedLocationCircles] = useState<any>(null);
  const [ionClass, setIonClass] = useState<string>('all');
  const [similarityScore, setSimilarityScore] = useState<number>(95);
  const [copiesNumber, setCopiesNumber] = useState<number>(1);
  const [showICG, setShowICG] = useState<boolean>(false);
  const [supermodel1, setSupermodel1] = useState<boolean>(true);
  const [supermodel2, setSupermodel2] = useState<boolean>(true);
  const [selectedNode, setSelectedNode] = useState<Node | null>(null);
  const [isLoading, setIsLoading] = useState<boolean>(true);
  const [isMobile, setIsMobile] = useState<boolean>(false);
  const [containerDimensions, setContainerDimensions] = useState<{width: number, height: number}>({width: 0, height: 0});
  
  // Check for mobile device and handle resize
  useEffect(() => {
    const checkMobile = () => {
      setIsMobile(window.innerWidth < 768); // Tailwind's md breakpoint
    };
    
    const handleResize = () => {
      checkMobile();
      // Trigger re-render of visualization by updating dependencies
      if (svgRef.current?.parentElement) {
        const container = svgRef.current.parentElement;
        const rect = container.getBoundingClientRect();
        setContainerDimensions({
          width: rect.width,
          height: rect.height
        });
      }
    };
    
    checkMobile();
    handleResize();
    window.addEventListener('resize', handleResize);
    return () => window.removeEventListener('resize', handleResize);
  }, []);

  // Load data on mount
  useEffect(() => {
    async function loadData() {
      try {
        // Load actual network data and fixed location circles
        const [networkResponse, fixedLocationResponse] = await Promise.all([
          fetch('/data/network_data.json'),
          fetch('/data/fixed_location_circles.json')
        ]);
        
        if (networkResponse.ok && fixedLocationResponse.ok) {
          const data = await networkResponse.json();
          const fixedLocations = await fixedLocationResponse.json();
          setFixedLocationCircles(fixedLocations);
          
          // Process the data to ensure it matches our interface
          const processedData: NetworkData = {
            nodes: data.nodes.map((node: any) => ({
              ...node,
              ion_class: node.ion_class || node.original_model?.ion_class || node.ion_channel_class || 'Other',
              icg: node.icg || node.original_model?.ICG || false,
              supermodel: node.supermodel || undefined,
              num_of_identicals: node.num_of_identicals || 1,
            })),
            links: data.links.map((link: any) => ({
              ...link,
              weight: link.weight || link.value || 95
            }))
          };
          
          setNetworkData(processedData);
        } else {
          console.error('Failed to load network data, using sample data');
          // Fallback to sample data
          const sampleData: NetworkData = {
            nodes: Array.from({ length: 100 }, (_, i) => ({
              id: `node${i}`,
              name: `Model ${i}`,
              ion_class: ['K', 'Na', 'Ca', 'IH', 'KCa', 'Other'][Math.floor(Math.random() * 6)],
              supermodel: Math.random() > 0.7 ? (Math.random() > 0.5 ? 1 : 2) : undefined,
              icg: Math.random() > 0.5,
            })),
            links: []
          };

          // Generate random links
          for (let i = 0; i < 150; i++) {
            const source = Math.floor(Math.random() * 100);
            const target = Math.floor(Math.random() * 100);
            if (source !== target) {
              sampleData.links.push({
                source: `node${source}`,
                target: `node${target}`,
                value: 75 + Math.random() * 25
              });
            }
          }

          setNetworkData(sampleData);
        }
        
        setIsLoading(false);
      } catch (error) {
        console.error('Error loading data:', error);
        setIsLoading(false);
      }
    }
    loadData();
  }, []);

  // D3 Force Simulation
  useEffect(() => {
    if (!networkData || !svgRef.current) return;

    // Get the actual container dimensions
    const container = svgRef.current.parentElement;
    if (!container) return;
    
    const containerRect = container.getBoundingClientRect();
    const width = containerRect.width - 32; // Account for padding
    const height = containerRect.height - 32; // Account for padding

    // Clear previous visualization
    d3.select(svgRef.current).selectAll('*').remove();

    const svg = d3.select(svgRef.current)
      .attr('width', width)
      .attr('height', height)
      .attr('viewBox', `0 0 ${width} ${height}`);

    const g = svg.append('g');

    // Filter data based on current filters
    const filteredNodes = networkData.nodes.filter(node => {
      if (ionClass !== 'all' && node.ion_class !== ionClass) return false;
      if (!showICG && node.icg) return false;
      if (!supermodel1 && node.supermodel === 1) return false;
      if (!supermodel2 && node.supermodel === 2) return false;
      return true;
    });

    const nodeIds = new Set(filteredNodes.map(n => n.id));
    const filteredLinks = networkData.links.filter(link => {
      const sourceId = typeof link.source === 'object' ? link.source.id : link.source;
      const targetId = typeof link.target === 'object' ? link.target.id : link.target;
      const linkWeight = link.weight || link.value || 0;
      return nodeIds.has(sourceId) && nodeIds.has(targetId) && linkWeight >= similarityScore;
    });

    // Group nodes based on filters
    const uniqueFilterCombinations = getUniqueFilterCombinations(ionClass, showICG, supermodel1, supermodel2);
    const nodesGroupedByFilter = assignNodesToGroups(filteredNodes, uniqueFilterCombinations, ionClass, showICG, supermodel1, supermodel2);
    const nodeIdToLocation = sortAndPositionGroups(nodesGroupedByFilter, fixedLocationCircles);
    
    // Check if we should split groups (more than one group)
    const splitVar = Object.keys(nodesGroupedByFilter).length > 1;

    // Force simulation constants (from original)
    const LINK_DISTANCE = 12;
    const CHARGE_STRENGTH_MULTIPLIER = -0.5;
    const CHARGE_STRENGTH_CONSTANT = -10.575;
    
    // Setup forces matching original implementation
    const linkForce = d3.forceLink(filteredLinks)
      .distance(() => LINK_DISTANCE)
      .id((d: any) => d.id);
    
    const chargeForce = d3.forceManyBody()
      .strength((d: any) => {
        const numIdenticals = d.num_of_identicals || 1;
        return (CHARGE_STRENGTH_MULTIPLIER * Math.pow(numIdenticals, 1.125) + CHARGE_STRENGTH_CONSTANT) + 
               (CHARGE_STRENGTH_MULTIPLIER * Math.pow(numIdenticals, 1.1275) + CHARGE_STRENGTH_CONSTANT);
      });
    
    const centerForce = d3.forceCenter(width / 2, height / 2);
    
    // Setup X and Y forces based on whether we're splitting groups
    const forceX = splitVar ? 
      d3.forceX((d: any) => nodeIdToLocation[d.id] ? nodeIdToLocation[d.id][0] * width : width / 2).strength(FORCE_STRENGTH) :
      d3.forceX(width / 2);
    
    const forceY = splitVar ?
      d3.forceY((d: any) => nodeIdToLocation[d.id] ? nodeIdToLocation[d.id][1] * height : height / 2).strength(FORCE_STRENGTH) :
      d3.forceY(height / 2);
    
    // Create simulation with original force configuration
    const simulation = d3.forceSimulation(filteredNodes)
      .force('links', linkForce)
      .force('charge', chargeForce)
      .force('center', centerForce)
      .force('x', forceX)
      .force('y', forceY);

    // Create links (using CSS for styling, matching original)
    const link = g.append('g')
      .attr('class', 'links')
      .selectAll('line')
      .data(filteredLinks)
      .enter().append('line')
      .attr('class', 'graph-link')
      .attr('stroke', '#aaa')
      .attr('stroke-opacity', 0.8);

    // Create nodes with dynamic radius based on num_of_identicals
    const node = g.append('g')
      .attr('class', 'nodes')
      .selectAll('circle')
      .data(filteredNodes)
      .enter().append('circle')
      .attr('r', (d: Node) => calculateNodeRadius(d))
      .attr('class', 'graph-node')
      .attr('fill', '#00BFFF')  // DeepSkyBlue - matching original
      .attr('stroke', '#aaa')
      .on('click', function(event: any, d: any) {
        // Reset all nodes to default color
        node.attr('fill', '#00BFFF')
            .attr('stroke', '#aaa')
            .classed('selected-node', false);
        
        // Highlight clicked node (matching original CSS classes)
        d3.select(this)
          .attr('fill', '#ffd700')  // Gold color for selected
          .attr('stroke', '#215885')
          .attr('stroke-width', 2)
          .classed('selected-node', true);
        
        setSelectedNode(d);
      })
      .call(d3.drag<any, any>()
        .on('start', dragstarted)
        .on('drag', dragged)
        .on('end', dragended) as any);

    // Add tooltips
    node.append('title')
      .text((d: Node) => `${d.name}\nClass: ${d.ion_class}\nICG: ${d.icg ? 'Yes' : 'No'}`);

    // Zoom behavior
    const zoom = d3.zoom()
      .scaleExtent([0.1, 10])
      .on('zoom', (event) => {
        g.attr('transform', event.transform);
      });

    svg.call(zoom as any);

    // Simulation tick (matching original - no boundary constraints)
    simulation.on('tick', () => {
      link
        .attr('x1', (d: any) => d.source.x)
        .attr('y1', (d: any) => d.source.y)
        .attr('x2', (d: any) => d.target.x)
        .attr('y2', (d: any) => d.target.y);

      node
        .attr('cx', (d: any) => d.x)
        .attr('cy', (d: any) => d.y);
    });

    // Drag functions
    function dragstarted(event: any, d: any) {
      if (!event.active) simulation.alphaTarget(0.3).restart();
      d.fx = d.x;
      d.fy = d.y;
    }

    function dragged(event: any, d: any) {
      d.fx = event.x;
      d.fy = event.y;
    }

    function dragended(event: any, d: any) {
      if (!event.active) simulation.alphaTarget(0);
      d.fx = null;
      d.fy = null;
    }

    return () => {
      simulation.stop();
    };
  }, [networkData, fixedLocationCircles, ionClass, similarityScore, showICG, supermodel1, supermodel2, copiesNumber, containerDimensions]);

  if (isLoading) {
    return (
      <div className="flex items-center justify-center h-screen bg-slate-50 dark:bg-slate-900">
        <div className="text-xl">Loading visualization...</div>
      </div>
    );
  }

  // Mobile warning
  if (isMobile) {
    return (
      <div className="flex items-center justify-center h-screen bg-slate-50 dark:bg-slate-900 p-8">
        <div className="text-center max-w-md">
          <svg className="w-24 h-24 mx-auto mb-6 text-slate-400" fill="none" stroke="currentColor" viewBox="0 0 24 24">
            <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={1.5} 
              d="M9.75 17L9 20l-1 1h8l-1-1-.75-3M3 13h18M5 17h14a2 2 0 002-2V5a2 2 0 00-2-2H5a2 2 0 00-2 2v10a2 2 0 002 2z" />
          </svg>
          <h2 className="text-2xl font-bold text-slate-900 dark:text-white mb-4">
            Desktop View Required
          </h2>
          <p className="text-slate-600 dark:text-slate-300 mb-6">
            The Ion Channel Network Visualizer requires a larger screen for the best experience. 
            Please access this tool on a tablet, laptop, or desktop computer.
          </p>
          <a
            href="/"
            className="inline-block px-6 py-3 rounded-lg bg-blue-600 text-white hover:bg-blue-700 transition-colors"
          >
            Return to Home
          </a>
        </div>
      </div>
    );
  }

  return (
    <div className="h-screen bg-slate-50 dark:bg-slate-900 overflow-hidden">
      {/* Header */}
      <div className="border-b border-slate-200 dark:border-slate-700 bg-white dark:bg-slate-800">
        <div className="px-4 py-3">
          <div className="flex justify-between items-center">
            <h1 className="text-xl font-bold text-slate-900 dark:text-white">
              Ion Channel Network Visualizer
            </h1>
            <a
              href="/"
              className="px-4 py-2 text-sm rounded-lg border border-slate-300 dark:border-slate-600 hover:bg-slate-100 dark:hover:bg-slate-700 transition-colors"
            >
              Back to Home
            </a>
          </div>
        </div>
      </div>

      <div className="flex h-[calc(100vh-57px)]">
        {/* Control Panel */}
        <div className="w-80 bg-white dark:bg-slate-800 border-r border-slate-200 dark:border-slate-700 overflow-y-auto">
          <div className="p-4 space-y-6">
            {/* Ion Channel Class Filter */}
            <div>
              <h3 className="text-sm font-semibold text-slate-700 dark:text-slate-300 mb-3">
                Ion Channel Class
              </h3>
              <div className="grid grid-cols-3 gap-2">
                {['all', 'K', 'Na', 'Ca', 'IH', 'KCa', 'Other'].map((cls) => (
                  <button
                    key={cls}
                    onClick={() => setIonClass(cls)}
                    className={`px-3 py-1.5 text-sm rounded-md transition-colors ${
                      ionClass === cls
                        ? 'bg-blue-600 text-white'
                        : 'bg-slate-100 dark:bg-slate-700 text-slate-700 dark:text-slate-300 hover:bg-slate-200 dark:hover:bg-slate-600'
                    }`}
                  >
                    {cls === 'all' ? 'All' : cls}
                  </button>
                ))}
              </div>
            </div>

            {/* Similarity Score */}
            <div>
              <h3 className="text-sm font-semibold text-slate-700 dark:text-slate-300 mb-3">
                Similarity Score: {similarityScore}%
              </h3>
              <div className="flex items-center space-x-3">
                <button
                  onClick={() => setSimilarityScore(Math.max(0, similarityScore - 5))}
                  className="w-8 h-8 rounded-md bg-slate-100 dark:bg-slate-700 hover:bg-slate-200 dark:hover:bg-slate-600 flex items-center justify-center"
                >
                  −
                </button>
                <input
                  type="range"
                  min="0"
                  max="100"
                  value={similarityScore}
                  onChange={(e) => setSimilarityScore(Number(e.target.value))}
                  className="flex-1"
                />
                <button
                  onClick={() => setSimilarityScore(Math.min(100, similarityScore + 5))}
                  className="w-8 h-8 rounded-md bg-slate-100 dark:bg-slate-700 hover:bg-slate-200 dark:hover:bg-slate-600 flex items-center justify-center"
                >
                  +
                </button>
              </div>
            </div>

            {/* Supermodels */}
            <div>
              <h3 className="text-sm font-semibold text-slate-700 dark:text-slate-300 mb-3">
                Supermodels
              </h3>
              <div className="space-y-2">
                <label className="flex items-center space-x-2">
                  <input
                    type="checkbox"
                    checked={supermodel1}
                    onChange={(e) => setSupermodel1(e.target.checked)}
                    className="rounded"
                  />
                  <span className="text-sm text-slate-600 dark:text-slate-400">Supermodel 1</span>
                </label>
                <label className="flex items-center space-x-2">
                  <input
                    type="checkbox"
                    checked={supermodel2}
                    onChange={(e) => setSupermodel2(e.target.checked)}
                    className="rounded"
                  />
                  <span className="text-sm text-slate-600 dark:text-slate-400">Supermodel 2</span>
                </label>
              </div>
            </div>

            {/* ICG Toggle */}
            <div>
              <label className="flex items-center justify-between">
                <span className="text-sm font-semibold text-slate-700 dark:text-slate-300">
                  Show ICG Entries
                </span>
                <button
                  onClick={() => setShowICG(!showICG)}
                  className={`relative inline-flex h-6 w-11 items-center rounded-full transition-colors ${
                    showICG ? 'bg-blue-600' : 'bg-gray-300 dark:bg-gray-600'
                  }`}
                >
                  <span
                    className={`inline-block h-4 w-4 transform rounded-full bg-white transition-transform ${
                      showICG ? 'translate-x-6' : 'translate-x-1'
                    }`}
                  />
                </button>
              </label>
            </div>

            {/* Selected Node Details */}
            {selectedNode && (
              <div className="border-t border-slate-200 dark:border-slate-700 pt-4">
                <h3 className="text-sm font-semibold text-slate-700 dark:text-slate-300 mb-3">
                  Selected Node
                </h3>
                <div className="bg-slate-50 dark:bg-slate-900 rounded-lg p-3 text-sm">
                  <p className="text-slate-600 dark:text-slate-400">
                    <span className="font-medium">Name:</span> {selectedNode.name}
                  </p>
                  <p className="text-slate-600 dark:text-slate-400">
                    <span className="font-medium">Class:</span> {selectedNode.ion_class}
                  </p>
                  <p className="text-slate-600 dark:text-slate-400">
                    <span className="font-medium">ICG:</span> {selectedNode.icg ? 'Yes' : 'No'}
                  </p>
                  {selectedNode.supermodel && (
                    <p className="text-slate-600 dark:text-slate-400">
                      <span className="font-medium">Supermodel:</span> {selectedNode.supermodel}
                    </p>
                  )}
                </div>
              </div>
            )}

            {/* Legend */}
            <div className="border-t border-slate-200 dark:border-slate-700 pt-4">
              <h3 className="text-sm font-semibold text-slate-700 dark:text-slate-300 mb-3">
                Legend
              </h3>
              <div className="space-y-2 text-sm">
                <div className="flex items-center space-x-2">
                  <div className="w-4 h-4 rounded-full" style={{ backgroundColor: '#00BFFF', border: '1px solid #aaa' }}></div>
                  <span className="text-slate-600 dark:text-slate-400">Ion Channel Model</span>
                </div>
                <div className="flex items-center space-x-2">
                  <div className="w-4 h-4 rounded-full" style={{ backgroundColor: '#ffd700', border: '1px solid #215885' }}></div>
                  <span className="text-slate-600 dark:text-slate-400">Selected/Source Node</span>
                </div>
                <div className="pt-2 text-xs text-slate-500 dark:text-slate-500">
                  Node size represents number of identical models
                </div>
              </div>
            </div>
          </div>
        </div>

        {/* Visualization Area */}
        <div className="flex-1 flex flex-col min-w-0">
          <div className="flex-1 bg-white dark:bg-slate-800 m-4 rounded-lg shadow-lg p-4 min-h-0">
            <svg 
              ref={svgRef} 
              className="w-full h-full border-2 border-slate-200 dark:border-slate-600 rounded"
              style={{ background: 'linear-gradient(to br, #fafafa, #f3f4f6)' }}
            ></svg>
          </div>
        </div>
      </div>
    </div>
  );
}
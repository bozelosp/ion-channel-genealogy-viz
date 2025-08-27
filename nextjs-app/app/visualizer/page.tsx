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

// Source code fetching function
async function fetchSourceCode(nodeIds: string[], networkData: NetworkData, fetchedFiles: {[key: string]: {source_code: string}}): Promise<{[key: string]: {source_code: string}}> {
  const baseUrl = "http://ion-channels.s3-website-eu-west-1.amazonaws.com/static/modelDB/";
  const updatedFetchedFiles = { ...fetchedFiles };

  for (const nodeId of nodeIds) {
    // Find the node data corresponding to the current node ID
    const nodeData = networkData.nodes.find(node => node.id === nodeId);
    if (!nodeData) {
      console.warn(`Node with ID ${nodeId} not found.`);
      continue;
    }

    const uniqueId = (nodeData as any).original_model?.unique_modelDB_mod_id;
    if (!uniqueId) {
      console.warn(`No unique_modelDB_mod_id found for node ${nodeId}`);
      continue;
    }

    // Check if the files for this node have already been fetched
    if (!updatedFetchedFiles[nodeId]) {
      const sourceCodePath = `${baseUrl}source_code/${uniqueId}.mod`;

      try {
        const response = await fetch(sourceCodePath);
        const sourceCode = await response.text();
        updatedFetchedFiles[nodeId] = { source_code: sourceCode };
      } catch (error) {
        console.error(`Error fetching file for node ${nodeId}:`, error);
      }
    }
  }

  return updatedFetchedFiles;
}

// Generate diff HTML from local API
async function generateDiffHtml(sourceCode: string, targetCode: string): Promise<string> {
  try {
    const diffUrl = `/api/generate-diff?string1=${encodeURIComponent(sourceCode)}&string2=${encodeURIComponent(targetCode)}`;
    const response = await fetch(diffUrl);
    if (!response.ok) {
      throw new Error(`HTTP error! status: ${response.status}`);
    }
    const diffHtml = await response.text();
    return diffHtml;
  } catch (error) {
    console.error('Diff generation error:', error);
    return '<div class="p-4 bg-red-50 border border-red-200 rounded-lg"><p class="text-red-700">Error generating diff: ' + (error instanceof Error ? error.message : 'Unknown error') + '</p></div>';
  }
}

// Handle diff creation and display
async function handleDiff(sourceIds: string[], targetIds: string[], fetchedFiles: {[key: string]: {source_code: string}}): Promise<string[]> {
  const diffBoxes: string[] = [];

  for (const sourceId of sourceIds) {
    const sourceCode = fetchedFiles[sourceId]?.source_code;
    if (!sourceCode) {
      console.warn(`Source code for node ${sourceId} is not fetched.`);
      continue;
    }

    for (const targetId of targetIds) {
      const targetCode = fetchedFiles[targetId]?.source_code;
      if (!targetCode) {
        console.warn(`Target code for node ${targetId} is not fetched.`);
        continue;
      }

      const diffHtml = await generateDiffHtml(sourceCode, targetCode);
      diffBoxes.push(diffHtml);
    }
  }

  return diffBoxes;
}

// Extend subgraph by connected nodes
function extendSubgraphByConnectedNodes(linkData: Link[], subgraphNodeIds: string[]): string[] {
  const extended = [...subgraphNodeIds];
  
  linkData.forEach(link => {
    const sourceId = typeof link.source === 'object' ? link.source.id : link.source;
    const targetId = typeof link.target === 'object' ? link.target.id : link.target;
    
    if (subgraphNodeIds.includes(sourceId) || subgraphNodeIds.includes(targetId)) {
      if (!extended.includes(sourceId)) extended.push(sourceId);
      if (!extended.includes(targetId)) extended.push(targetId);
    }
  });
  
  return [...new Set(extended)]; // Remove duplicates
}

// Generate all combinations between source and target nodes
function generateAllCombinations(sourceNodes: Node[], targetNodes: Node[]): {source: Node, target: Node}[] {
  const combinations: {source: Node, target: Node}[] = [];
  
  for (const source of sourceNodes) {
    for (const target of targetNodes) {
      if (source.id !== target.id) {
        combinations.push({ source, target });
      }
    }
  }
  
  return combinations;
}

// Get selected nodes by class
function getSelectedNodesByClass(svg: any): {sources: Node[], selected: Node[], all: Node[]} {
  const sources: Node[] = [];
  const selected: Node[] = [];
  const all: Node[] = [];
  
  svg.selectAll('.source-node').each((d: Node) => {
    sources.push(d);
    all.push(d);
  });
  
  svg.selectAll('.selected-node').each((d: Node) => {
    if (!sources.find(s => s.id === d.id)) {
      selected.push(d);
      all.push(d);
    }
  });
  
  return { sources, selected, all };
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
  
  // Source code comparison states
  const [sourceNodeIds, setSourceNodeIds] = useState<string[]>([]);
  const [targetNodeIds, setTargetNodeIds] = useState<string[]>([]);
  const [fetchedFiles, setFetchedFiles] = useState<{[key: string]: {source_code: string}}>({});
  const [diffBoxes, setDiffBoxes] = useState<string[]>([]);
  const [currentDiffIndex, setCurrentDiffIndex] = useState<number>(0);
  const [showDiffView, setShowDiffView] = useState<boolean>(false);
  
  // View control
  const [fitToViewFunction, setFitToViewFunction] = useState<(() => void) | null>(null);
  
  // Context menu state
  const [contextMenu, setContextMenu] = useState<{x: number, y: number, node: Node} | null>(null);
  const [diffCombinations, setDiffCombinations] = useState<{source: Node, target: Node, html?: string}[]>([]);
  const [currentCombinationIndex, setCurrentCombinationIndex] = useState<number>(0);
  const [isGeneratingDiffs, setIsGeneratingDiffs] = useState<boolean>(false);
  
  // Group summary state
  const [groupSummaries, setGroupSummaries] = useState<{key: string, nodeCount: number}[]>([]);
  
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

  // Keyboard event handler for shortcuts
  useEffect(() => {
    const handleKeyPress = async (event: KeyboardEvent) => {
      if (event.code === 'KeyD') {
        // Generate awesome diffs when 'D' key is pressed
        if (svgRef.current) {
          const svg = d3.select(svgRef.current);
          const nodesByClass = getSelectedNodesByClass(svg);
          const allSelectedNodes = nodesByClass.all;
          
          if (allSelectedNodes.length >= 2) {
            await generateAwesomeDiffs(allSelectedNodes);
          } else {
            alert('Please select at least 2 nodes to generate diffs');
          }
        }
      } else if (event.code === 'ArrowLeft' && showDiffView && diffCombinations.length > 0) {
        // Navigate to previous diff combination
        setCurrentCombinationIndex(prev => Math.max(0, prev - 1));
      } else if (event.code === 'ArrowRight' && showDiffView && diffCombinations.length > 0) {
        // Navigate to next diff combination
        setCurrentCombinationIndex(prev => Math.min(diffCombinations.length - 1, prev + 1));
      } else if (event.code === 'KeyF') {
        // Fit to view when 'F' key is pressed
        if (fitToViewFunction) {
          fitToViewFunction();
        }
      } else if (event.code === 'Escape') {
        // Close diff view and context menu
        setShowDiffView(false);
        setDiffCombinations([]);
        setCurrentCombinationIndex(0);
        setContextMenu(null);
      }
    };

    window.addEventListener('keydown', handleKeyPress);
    return () => window.removeEventListener('keydown', handleKeyPress);
  }, [showDiffView, diffCombinations]);

  // Close context menu when clicking elsewhere
  useEffect(() => {
    const handleClickOutside = () => setContextMenu(null);
    
    if (contextMenu) {
      document.addEventListener('click', handleClickOutside);
      return () => document.removeEventListener('click', handleClickOutside);
    }
  }, [contextMenu]);

  // Awesome diff generation function
  const generateAwesomeDiffs = async (selectedNodes: Node[]) => {
    if (!networkData) return;
    
    setIsGeneratingDiffs(true);
    
    // Determine source and target nodes based on selection
    let sourceNodes: Node[] = [];
    let targetNodes: Node[] = [];
    
    if (svgRef.current) {
      const svg = d3.select(svgRef.current);
      const nodesByClass = getSelectedNodesByClass(svg);
      
      if (nodesByClass.sources.length > 0) {
        // Use designated source nodes
        sourceNodes = nodesByClass.sources;
        targetNodes = nodesByClass.selected.length > 0 ? nodesByClass.selected : selectedNodes.filter(n => !sourceNodes.find(s => s.id === n.id));
      } else {
        // No designated sources, use all combinations
        sourceNodes = selectedNodes;
        targetNodes = selectedNodes;
      }
    }
    
    // Generate all combinations
    const combinations = generateAllCombinations(sourceNodes, targetNodes);
    
    if (combinations.length === 0) {
      alert('Please select at least 2 different nodes or designate source nodes with Ctrl+Click');
      setIsGeneratingDiffs(false);
      return;
    }
    
    // Fetch source code for all involved nodes
    const allNodeIds = [...new Set([...sourceNodes.map(n => n.id), ...targetNodes.map(n => n.id)])];
    const updatedFiles = await fetchSourceCode(allNodeIds, networkData, fetchedFiles);
    setFetchedFiles(updatedFiles);
    
    // Generate diffs for each combination
    const combinationsWithDiffs: {source: Node, target: Node, html?: string}[] = [];
    
    for (const combo of combinations) {
      const sourceCode = updatedFiles[combo.source.id]?.source_code;
      const targetCode = updatedFiles[combo.target.id]?.source_code;
      
      if (sourceCode && targetCode) {
        const diffHtml = await generateDiffHtml(sourceCode, targetCode);
        combinationsWithDiffs.push({
          ...combo,
          html: diffHtml
        });
      } else {
        combinationsWithDiffs.push(combo);
      }
    }
    
    setDiffCombinations(combinationsWithDiffs);
    setCurrentCombinationIndex(0);
    setShowDiffView(true);
    setIsGeneratingDiffs(false);
    setContextMenu(null);
  };

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
    
    // Create group summaries (matching original implementation)
    const newGroupSummaries: {key: string, nodeCount: number}[] = [];
    Object.entries(nodesGroupedByFilter).forEach(([groupKey, groupNodes]) => {
      if (groupNodes.length > 0) {
        // Format the group key like the original (replace commas, format nicely)
        let formattedKey = groupKey
          .replace(',', ': ')
          .replace(/,/g, ' • ')
          .replace('Supermodel 1 • Supermodel 2', 'Supermodel 1 & 2');
        
        newGroupSummaries.push({
          key: formattedKey,
          nodeCount: groupNodes.length
        });
      }
    });
    setGroupSummaries(newGroupSummaries);
    
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
      .on('click', async function(event: any, d: any) {
        const currentNode = d3.select(this);

        if (event.shiftKey) {
          // Shift+Click: Select subgraph and find source node
          let subgraphNodeIds = [d.id];
          
          let previousSize = -1;
          let currentSize = 1;
          
          // Expand subgraph until no new nodes are added
          while (previousSize !== currentSize) {
            previousSize = currentSize;
            subgraphNodeIds = extendSubgraphByConnectedNodes(filteredLinks, subgraphNodeIds);
            currentSize = subgraphNodeIds.length;
          }
          
          // Reset all node classes
          node.classed('selected-node', false).classed('source-node', false)
              .attr('fill', '#00BFFF').attr('stroke', '#aaa').attr('stroke-width', 1);
          
          // Get nodes in subgraph
          const nodesInSubgraph = node.filter((nodeData: any) => subgraphNodeIds.includes(nodeData.id));
          nodesInSubgraph.classed('selected-node', true)
                         .attr('stroke', '#215885')
                         .attr('stroke-width', 2);
          
          // Find the earliest year among nodes in subgraph
          const subgraphData = nodesInSubgraph.data();
          const yearsInSubgraph = subgraphData.map((nodeData: any) => 
            nodeData.original_model?.Year || 2023
          );
          const earliestYear = Math.min(...yearsInSubgraph);
          
          // Filter nodes by earliest year
          let originalNodes = subgraphData.filter((nodeData: any) => 
            (nodeData.original_model?.Year || 2023) === earliestYear
          );
          
          // If multiple nodes have same earliest year, choose by smallest modelDB_dir
          if (originalNodes.length > 1) {
            const modelDBDirs = originalNodes.map((nodeData: any) => 
              parseInt(nodeData.original_model?.modelDB_dir || '999999')
            );
            const smallestModelDBDir = Math.min(...modelDBDirs);
            originalNodes = originalNodes.filter((nodeData: any) => 
              parseInt(nodeData.original_model?.modelDB_dir || '999999') === smallestModelDBDir
            );
          }
          
          // Mark source node(s)
          const sourceIds = originalNodes.map((nodeData: any) => nodeData.id);
          node.filter((nodeData: any) => sourceIds.includes(nodeData.id))
              .classed('source-node', true)
              .attr('fill', '#ffd700')
              .attr('stroke', '#215885')
              .attr('stroke-width', 2);
          
          // Set state for diff generation
          const newSourceIds = sourceIds;
          const newTargetIds = subgraphNodeIds.filter(id => !sourceIds.includes(id));
          
          setSourceNodeIds(newSourceIds);
          setTargetNodeIds(newTargetIds);
          
          // Fetch source code for all nodes in subgraph
          if (networkData) {
            const allNodeIds = [...newSourceIds, ...newTargetIds];
            const updatedFetchedFiles = await fetchSourceCode(allNodeIds, networkData, fetchedFiles);
            setFetchedFiles(updatedFetchedFiles);
          }
          
        } else if (event.ctrlKey || event.metaKey) {
          // Ctrl/Cmd+Click: Toggle source node
          const isSourceNode = currentNode.classed('source-node');
          
          if (isSourceNode) {
            currentNode.classed('source-node', false)
                       .attr('fill', '#00BFFF')
                       .attr('stroke', '#aaa')
                       .attr('stroke-width', 1);
          } else {
            currentNode.classed('source-node', true)
                       .attr('fill', '#ffd700')
                       .attr('stroke', '#215885')
                       .attr('stroke-width', 2);
          }
          
        } else {
          // Regular click: Toggle selected node
          const isSelectedNode = currentNode.classed('selected-node');
          
          if (isSelectedNode) {
            currentNode.classed('selected-node', false)
                       .attr('fill', '#00BFFF')
                       .attr('stroke', '#aaa')
                       .attr('stroke-width', 1);
            setSelectedNode(null);
          } else {
            currentNode.classed('selected-node', true)
                       .attr('fill', '#00BFFF')
                       .attr('stroke', '#215885')
                       .attr('stroke-width', 2);
            setSelectedNode(d);
          }
        }
      })
      .on('contextmenu', function(event: any, d: any) {
        event.preventDefault();
        
        // Get mouse position relative to the page
        const rect = (event.target as Element).getBoundingClientRect();
        const x = event.clientX;
        const y = event.clientY;
        
        setContextMenu({ x, y, node: d });
      })
      .call(d3.drag<any, any>()
        .on('start', dragstarted)
        .on('drag', dragged)
        .on('end', dragended) as any);

    // Add tooltips
    node.append('title')
      .text((d: Node) => `${d.name}\nClass: ${d.ion_class}\nICG: ${d.icg ? 'Yes' : 'No'}`);

    // Calculate bounding box of all nodes for proper zoom setup
    const calculateBounds = () => {
      let minX = Infinity, maxX = -Infinity;
      let minY = Infinity, maxY = -Infinity;
      let hasValidBounds = false;
      
      filteredNodes.forEach(node => {
        const nodeRadius = calculateNodeRadius(node);
        // Use actual node positions if available, otherwise fall back to center
        const x = (typeof node.x === 'number' && !isNaN(node.x)) ? node.x : width / 2;
        const y = (typeof node.y === 'number' && !isNaN(node.y)) ? node.y : height / 2;
        
        minX = Math.min(minX, x - nodeRadius);
        maxX = Math.max(maxX, x + nodeRadius);
        minY = Math.min(minY, y - nodeRadius);
        maxY = Math.max(maxY, y + nodeRadius);
        hasValidBounds = true;
      });
      
      // If we don't have valid bounds, use reasonable defaults
      if (!hasValidBounds || minX === Infinity) {
        const defaultRadius = 100;
        minX = width / 2 - defaultRadius;
        maxX = width / 2 + defaultRadius;
        minY = height / 2 - defaultRadius;
        maxY = height / 2 + defaultRadius;
      }
      
      return { minX, maxX, minY, maxY };
    };

    // Setup zoom behavior with proper bounds and centering
    const zoom = d3.zoom()
      .scaleExtent([0.05, 20])  // Allow more zoom out and zoom in
      .on('zoom', (event) => {
        g.attr('transform', event.transform);
      });

    svg.call(zoom as any);

    // Function to fit visualization with aesthetic margins
    const fitToView = () => {
      const bounds = calculateBounds();
      const boundsWidth = bounds.maxX - bounds.minX;
      const boundsHeight = bounds.maxY - bounds.minY;
      
      // Add aesthetic margins (20% of container size)
      const margin = 0.2;
      const marginX = width * margin;
      const marginY = height * margin;
      
      // Available space for content
      const availableWidth = width - (2 * marginX);
      const availableHeight = height - (2 * marginY);
      
      // Calculate scale to fit content with margins
      const scale = Math.min(
        availableWidth / boundsWidth,
        availableHeight / boundsHeight,
        1 // Don't zoom in beyond 1:1 initially
      );
      
      // Center point of the bounds
      const centerX = (bounds.minX + bounds.maxX) / 2;
      const centerY = (bounds.minY + bounds.maxY) / 2;
      
      // Calculate translation to center the content
      const translateX = width / 2 - centerX * scale;
      const translateY = height / 2 - centerY * scale;
      
      // Apply the transform smoothly
      const transform = d3.zoomIdentity
        .translate(translateX, translateY)
        .scale(scale);
        
      svg.transition()
        .duration(750)
        .call(zoom.transform, transform);
    };

    // Store the fit function for keyboard access
    setFitToViewFunction(() => fitToView);

    // Initial fit after simulation settles (adaptive timing)
    let fittingTimeout: NodeJS.Timeout;
    const scheduleInitialFit = () => {
      clearTimeout(fittingTimeout);
      fittingTimeout = setTimeout(() => {
        fitToView();
      }, 1500);
    };
    
    // Schedule initial fit and re-schedule if simulation is still active
    scheduleInitialFit();
    
    // Re-schedule fit if simulation restarts (e.g., due to interactions)
    simulation.on('end', () => {
      setTimeout(fitToView, 300);
    });

    // Add fit-to-view functionality on double-click
    svg.on('dblclick.zoom', null); // Remove default double-click zoom
    svg.on('dblclick', fitToView);

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
              className="px-4 py-2 text-sm rounded-lg border border-slate-300 dark:border-slate-600 hover:bg-slate-100 dark:hover:bg-slate-700 transition-colors cursor-pointer"
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
                    className={`px-3 py-1.5 text-sm rounded-md transition-colors cursor-pointer ${
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
                  className="w-8 h-8 rounded-md bg-slate-100 dark:bg-slate-700 hover:bg-slate-200 dark:hover:bg-slate-600 flex items-center justify-center cursor-pointer"
                >
                  −
                </button>
                <input
                  type="range"
                  min="0"
                  max="100"
                  value={similarityScore}
                  onChange={(e) => setSimilarityScore(Number(e.target.value))}
                  className="flex-1 cursor-pointer"
                />
                <button
                  onClick={() => setSimilarityScore(Math.min(100, similarityScore + 5))}
                  className="w-8 h-8 rounded-md bg-slate-100 dark:bg-slate-700 hover:bg-slate-200 dark:hover:bg-slate-600 flex items-center justify-center cursor-pointer"
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
                <label className="flex items-center space-x-2 cursor-pointer">
                  <input
                    type="checkbox"
                    checked={supermodel1}
                    onChange={(e) => setSupermodel1(e.target.checked)}
                    className="rounded cursor-pointer"
                  />
                  <span className="text-sm text-slate-600 dark:text-slate-400">Supermodel 1</span>
                </label>
                <label className="flex items-center space-x-2 cursor-pointer">
                  <input
                    type="checkbox"
                    checked={supermodel2}
                    onChange={(e) => setSupermodel2(e.target.checked)}
                    className="rounded cursor-pointer"
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
                  className={`relative inline-flex h-6 w-11 items-center rounded-full transition-colors cursor-pointer ${
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

            {/* Keyboard Shortcuts */}
            <div className="border-t border-slate-200 dark:border-slate-700 pt-4">
              <h3 className="text-sm font-semibold text-slate-700 dark:text-slate-300 mb-3">
                Interaction Guide
              </h3>
              <div className="space-y-3 text-xs text-slate-600 dark:text-slate-400">
                <div className="space-y-1">
                  <div className="font-medium text-slate-700 dark:text-slate-300">Node Selection:</div>
                  <div><kbd className="px-1.5 py-0.5 bg-slate-200 dark:bg-slate-700 rounded text-xs">Click</kbd> Select/deselect</div>
                  <div><kbd className="px-1.5 py-0.5 bg-slate-200 dark:bg-slate-700 rounded text-xs">Shift+Click</kbd> Select subgraph</div>
                  <div><kbd className="px-1.5 py-0.5 bg-slate-200 dark:bg-slate-700 rounded text-xs">Ctrl+Click</kbd> Toggle source</div>
                  <div><kbd className="px-1.5 py-0.5 bg-slate-200 dark:bg-slate-700 rounded text-xs">Right-Click</kbd> Context menu</div>
                </div>
                <div className="space-y-1">
                  <div className="font-medium text-slate-700 dark:text-slate-300">Code Comparison:</div>
                  <div><kbd className="px-1.5 py-0.5 bg-slate-200 dark:bg-slate-700 rounded text-xs">D</kbd> Generate all combinations</div>
                  <div><kbd className="px-1.5 py-0.5 bg-slate-200 dark:bg-slate-700 rounded text-xs">←/→</kbd> Navigate combinations</div>
                  <div><kbd className="px-1.5 py-0.5 bg-slate-200 dark:bg-slate-700 rounded text-xs">Esc</kbd> Close diff view</div>
                </div>
                <div className="space-y-1">
                  <div className="font-medium text-slate-700 dark:text-slate-300">View Control:</div>
                  <div><kbd className="px-1.5 py-0.5 bg-slate-200 dark:bg-slate-700 rounded text-xs">F</kbd> Fit to view</div>
                  <div><kbd className="px-1.5 py-0.5 bg-slate-200 dark:bg-slate-700 rounded text-xs">Double-Click</kbd> Fit to view</div>
                </div>
              </div>
            </div>
          </div>
        </div>

        {/* Visualization Area */}
        <div className="flex-1 flex flex-col min-w-0 relative">
          <div className="flex-1 bg-white dark:bg-slate-800 m-4 rounded-lg shadow-lg p-4 min-h-0">
            <svg 
              ref={svgRef} 
              className="w-full h-full border-2 border-slate-200 dark:border-slate-600 rounded"
              style={{ background: 'linear-gradient(to br, #fafafa, #f3f4f6)' }}
            ></svg>
          </div>

          {/* Group Summary Boxes (matching original position and styling) */}
          {groupSummaries.length > 1 && (
            <div 
              className="absolute overflow-y-auto rounded-md"
              style={{
                right: '10px',
                top: '282px', 
                width: '354px',
                maxHeight: 'calc(100vh - 320px)'
              }}
            >
                {groupSummaries.map((group, index) => (
                  <div 
                    key={index}
                    className="group-summary-box"
                  >
                    <span className="group-summary-box-title">
                      {group.key}
                    </span>
                    <br />
                    <span className="group-summary-box-node-count">
                      {group.nodeCount} nodes
                    </span>
                  </div>
                ))}
            </div>
          )}
        </div>
      </div>

      {/* Context Menu */}
      {contextMenu && (
        <div 
          className="fixed z-50 bg-white dark:bg-slate-800 rounded-lg shadow-xl border border-slate-200 dark:border-slate-700 py-2 min-w-48"
          style={{
            left: contextMenu.x,
            top: contextMenu.y,
            transform: 'translate(-50%, -10px)'
          }}
          onClick={(e) => e.stopPropagation()}
        >
          <div className="px-4 py-2 border-b border-slate-200 dark:border-slate-700">
            <div className="text-sm font-semibold text-slate-900 dark:text-white">
              {contextMenu.node.name || contextMenu.node.id}
            </div>
            <div className="text-xs text-slate-500 dark:text-slate-400">
              {contextMenu.node.ion_class} • {contextMenu.node.icg ? 'ICG' : 'ModelDB'}
            </div>
          </div>
          
          <button
            className="w-full px-4 py-2 text-left text-sm text-slate-700 dark:text-slate-300 hover:bg-slate-100 dark:hover:bg-slate-700 flex items-center space-x-2 cursor-pointer"
            onClick={async () => {
              if (svgRef.current) {
                const svg = d3.select(svgRef.current);
                const nodesByClass = getSelectedNodesByClass(svg);
                const allSelectedNodes = [...nodesByClass.all, contextMenu.node];
                await generateAwesomeDiffs(allSelectedNodes);
              }
            }}
            disabled={isGeneratingDiffs}
          >
            <svg className="w-4 h-4" fill="none" stroke="currentColor" viewBox="0 0 24 24">
              <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2} d="M9 12h6m-6 4h6m2 5H7a2 2 0 01-2-2V5a2 2 0 012-2h5.586a1 1 0 01.707.293l5.414 5.414a1 1 0 01.293.707V19a2 2 0 01-2 2z" />
            </svg>
            <span>{isGeneratingDiffs ? 'Generating...' : 'Compare with selected nodes'}</span>
          </button>
          
          <button
            className="w-full px-4 py-2 text-left text-sm text-slate-700 dark:text-slate-300 hover:bg-slate-100 dark:hover:bg-slate-700 flex items-center space-x-2 cursor-pointer"
            onClick={() => {
              setContextMenu(null);
              // Toggle node selection
              if (svgRef.current) {
                const svg = d3.select(svgRef.current);
                const targetNode = svg.selectAll('circle').filter((d: any) => d.id === contextMenu.node.id);
                const isSelected = targetNode.classed('selected-node');
                
                if (isSelected) {
                  targetNode.classed('selected-node', false)
                             .attr('fill', '#00BFFF')
                             .attr('stroke', '#aaa')
                             .attr('stroke-width', 1);
                } else {
                  targetNode.classed('selected-node', true)
                             .attr('fill', '#00BFFF')
                             .attr('stroke', '#215885')
                             .attr('stroke-width', 2);
                }
              }
            }}
          >
            <svg className="w-4 h-4" fill="none" stroke="currentColor" viewBox="0 0 24 24">
              <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2} d="M5 13l4 4L19 7" />
            </svg>
            <span>Toggle selection</span>
          </button>

          <button
            className="w-full px-4 py-2 text-left text-sm text-slate-700 dark:text-slate-300 hover:bg-slate-100 dark:hover:bg-slate-700 flex items-center space-x-2 cursor-pointer"
            onClick={() => {
              setContextMenu(null);
              // Toggle source node
              if (svgRef.current) {
                const svg = d3.select(svgRef.current);
                const targetNode = svg.selectAll('circle').filter((d: any) => d.id === contextMenu.node.id);
                const isSource = targetNode.classed('source-node');
                
                if (isSource) {
                  targetNode.classed('source-node', false)
                             .attr('fill', '#00BFFF')
                             .attr('stroke', '#aaa')
                             .attr('stroke-width', 1);
                } else {
                  targetNode.classed('source-node', true)
                             .attr('fill', '#ffd700')
                             .attr('stroke', '#215885')
                             .attr('stroke-width', 2);
                }
              }
            }}
          >
            <svg className="w-4 h-4" fill="none" stroke="currentColor" viewBox="0 0 24 24">
              <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2} d="M11.049 2.927c.3-.921 1.603-.921 1.902 0l1.519 4.674a1 1 0 00.95.69h4.915c.969 0 1.371 1.24.588 1.81l-3.976 2.888a1 1 0 00-.363 1.118l1.518 4.674c.3.922-.755 1.688-1.538 1.118l-3.976-2.888a1 1 0 00-1.176 0l-3.976 2.888c-.783.57-1.838-.197-1.538-1.118l1.518-4.674a1 1 0 00-.363-1.118l-3.976-2.888c-.784-.57-.38-1.81.588-1.81h4.914a1 1 0 00.951-.69l1.519-4.674z" />
            </svg>
            <span>Toggle as source</span>
          </button>
        </div>
      )}

      {/* Awesome Diff View Overlay */}
      {showDiffView && diffCombinations.length > 0 && (
        <div className="fixed inset-0 bg-black bg-opacity-75 z-50 flex items-center justify-center p-4">
          <div className="bg-white dark:bg-slate-800 rounded-xl shadow-2xl max-w-7xl w-full max-h-[95vh] overflow-hidden border border-slate-200 dark:border-slate-700">
            {/* Awesome Diff Header */}
            <div className="bg-gradient-to-r from-blue-600 to-purple-600 text-white p-6">
              <div className="flex items-center justify-between">
                <div className="flex items-center space-x-4">
                  <div className="bg-white bg-opacity-20 rounded-lg p-3">
                    <svg className="w-6 h-6" fill="none" stroke="currentColor" viewBox="0 0 24 24">
                      <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2} d="M9 12h6m-6 4h6m2 5H7a2 2 0 01-2-2V5a2 2 0 012-2h5.586a1 1 0 01.707.293l5.414 5.414a1 1 0 01.293.707V19a2 2 0 01-2 2z" />
                    </svg>
                  </div>
                  <div>
                    <h3 className="text-2xl font-bold">Code Comparison Matrix</h3>
                    <p className="text-blue-100 text-sm">
                      Exploring {diffCombinations.length} unique combinations
                    </p>
                  </div>
                </div>
                <button
                  onClick={() => {
                    setShowDiffView(false);
                    setDiffCombinations([]);
                    setCurrentCombinationIndex(0);
                  }}
                  className="bg-white bg-opacity-20 hover:bg-opacity-30 rounded-lg p-2 transition-all duration-200 cursor-pointer"
                >
                  <svg className="w-6 h-6" fill="none" stroke="currentColor" viewBox="0 0 24 24">
                    <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2} d="M6 18L18 6M6 6l12 12" />
                  </svg>
                </button>
              </div>
            </div>

            {/* Current Comparison Info */}
            <div className="bg-slate-50 dark:bg-slate-900 border-b border-slate-200 dark:border-slate-700 p-4">
              <div className="flex items-center justify-between">
                <div className="flex items-center space-x-6">
                  <div className="flex items-center space-x-2">
                    <div className="w-3 h-3 rounded-full bg-gradient-to-r from-yellow-400 to-orange-500"></div>
                    <span className="text-sm font-medium text-slate-700 dark:text-slate-300">Source:</span>
                    <span className="text-sm text-slate-600 dark:text-slate-400 font-mono">
                      {diffCombinations[currentCombinationIndex]?.source.name || diffCombinations[currentCombinationIndex]?.source.id}
                    </span>
                    <span className="text-xs text-slate-500 dark:text-slate-500 bg-slate-200 dark:bg-slate-700 px-2 py-1 rounded">
                      {diffCombinations[currentCombinationIndex]?.source.ion_class}
                    </span>
                  </div>
                  <div className="text-slate-400">→</div>
                  <div className="flex items-center space-x-2">
                    <div className="w-3 h-3 rounded-full bg-gradient-to-r from-blue-400 to-cyan-500"></div>
                    <span className="text-sm font-medium text-slate-700 dark:text-slate-300">Target:</span>
                    <span className="text-sm text-slate-600 dark:text-slate-400 font-mono">
                      {diffCombinations[currentCombinationIndex]?.target.name || diffCombinations[currentCombinationIndex]?.target.id}
                    </span>
                    <span className="text-xs text-slate-500 dark:text-slate-500 bg-slate-200 dark:bg-slate-700 px-2 py-1 rounded">
                      {diffCombinations[currentCombinationIndex]?.target.ion_class}
                    </span>
                  </div>
                </div>
                
                <div className="flex items-center space-x-3">
                  <span className="text-sm text-slate-500 dark:text-slate-400">
                    {currentCombinationIndex + 1} of {diffCombinations.length}
                  </span>
                  <div className="flex items-center space-x-1">
                    <button
                      onClick={() => setCurrentCombinationIndex(Math.max(0, currentCombinationIndex - 1))}
                      disabled={currentCombinationIndex === 0}
                      className="p-2 rounded-lg bg-slate-200 dark:bg-slate-700 hover:bg-slate-300 dark:hover:bg-slate-600 disabled:opacity-50 disabled:cursor-not-allowed transition-all duration-200"
                    >
                      <svg className="w-4 h-4" fill="none" stroke="currentColor" viewBox="0 0 24 24">
                        <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2} d="M15 19l-7-7 7-7" />
                      </svg>
                    </button>
                    <button
                      onClick={() => setCurrentCombinationIndex(Math.min(diffCombinations.length - 1, currentCombinationIndex + 1))}
                      disabled={currentCombinationIndex === diffCombinations.length - 1}
                      className="p-2 rounded-lg bg-slate-200 dark:bg-slate-700 hover:bg-slate-300 dark:hover:bg-slate-600 disabled:opacity-50 disabled:cursor-not-allowed transition-all duration-200"
                    >
                      <svg className="w-4 h-4" fill="none" stroke="currentColor" viewBox="0 0 24 24">
                        <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2} d="M9 5l7 7-7 7" />
                      </svg>
                    </button>
                  </div>
                </div>
              </div>
            </div>

            {/* Progress Bar */}
            <div className="bg-slate-100 dark:bg-slate-800 h-1">
              <div 
                className="h-full bg-gradient-to-r from-blue-500 to-purple-500 transition-all duration-300 ease-out"
                style={{ width: `${((currentCombinationIndex + 1) / diffCombinations.length) * 100}%` }}
              />
            </div>

            {/* Diff Content */}
            <div className="flex-1 overflow-auto max-h-[calc(95vh-200px)]">
              {diffCombinations[currentCombinationIndex]?.html ? (
                <div className="p-6">
                  <div 
                    className="diff-container prose max-w-none"
                    dangerouslySetInnerHTML={{ __html: diffCombinations[currentCombinationIndex].html }}
                    style={{ 
                      fontFamily: 'ui-monospace, SFMono-Regular, "SF Mono", Monaco, Consolas, "Liberation Mono", "Courier New", monospace',
                      fontSize: '13px',
                      lineHeight: '1.6'
                    }}
                  />
                </div>
              ) : (
                <div className="flex items-center justify-center h-64">
                  <div className="text-center">
                    <div className="animate-spin rounded-full h-12 w-12 border-b-2 border-blue-600 mx-auto mb-4"></div>
                    <p className="text-slate-600 dark:text-slate-400">Loading comparison...</p>
                  </div>
                </div>
              )}
            </div>

            {/* Quick Navigation */}
            <div className="bg-slate-50 dark:bg-slate-900 border-t border-slate-200 dark:border-slate-700 p-4">
              <div className="flex items-center justify-between text-sm text-slate-600 dark:text-slate-400">
                <div className="flex items-center space-x-4">
                  <kbd className="px-2 py-1 bg-slate-200 dark:bg-slate-700 rounded text-xs">←</kbd>
                  <span>Previous</span>
                  <kbd className="px-2 py-1 bg-slate-200 dark:bg-slate-700 rounded text-xs">→</kbd>
                  <span>Next</span>
                  <kbd className="px-2 py-1 bg-slate-200 dark:bg-slate-700 rounded text-xs">Esc</kbd>
                  <span>Close</span>
                </div>
                
                <div className="flex items-center space-x-2">
                  <span className="text-xs">Jump to:</span>
                  <select
                    value={currentCombinationIndex}
                    onChange={(e) => setCurrentCombinationIndex(Number(e.target.value))}
                    className="text-xs bg-white dark:bg-slate-800 border border-slate-300 dark:border-slate-600 rounded px-2 py-1"
                  >
                    {diffCombinations.map((combo, index) => (
                      <option key={index} value={index}>
                        {combo.source.name || combo.source.id} → {combo.target.name || combo.target.id}
                      </option>
                    ))}
                  </select>
                </div>
              </div>
            </div>
          </div>
        </div>
      )}
    </div>
  );
}
'use client';

import { useEffect, useRef, useState } from 'react';
import * as d3 from 'd3';

interface Node {
  id: string;
  name: string;
  ion_class?: string;
  supermodel?: number;
  icg?: boolean;
  x?: number;
  y?: number;
  fx?: number | null;
  fy?: number | null;
}

interface Link {
  source: string | Node;
  target: string | Node;
  value: number;
}

interface NetworkData {
  nodes: Node[];
  links: Link[];
}

export default function Visualizer() {
  const svgRef = useRef<SVGSVGElement>(null);
  const [networkData, setNetworkData] = useState<NetworkData | null>(null);
  const [ionClass, setIonClass] = useState<string>('all');
  const [similarityScore, setSimilarityScore] = useState<number>(95);
  const [copiesNumber, setCopiesNumber] = useState<number>(1);
  const [showICG, setShowICG] = useState<boolean>(false);
  const [supermodel1, setSupermodel1] = useState<boolean>(true);
  const [supermodel2, setSupermodel2] = useState<boolean>(true);
  const [selectedNode, setSelectedNode] = useState<Node | null>(null);
  const [isLoading, setIsLoading] = useState<boolean>(true);

  // Load data on mount
  useEffect(() => {
    async function loadData() {
      try {
        // Load actual network data
        const response = await fetch('/data/network_data.json');
        if (response.ok) {
          const data = await response.json();
          
          // Process the data to ensure it matches our interface
          const processedData: NetworkData = {
            nodes: data.nodes.map((node: any) => ({
              ...node,
              ion_class: node.ion_class || node.ion_channel_class || 'Other',
              icg: node.icg || false,
              supermodel: node.supermodel || undefined,
            })),
            links: data.links.map((link: any) => ({
              ...link,
              value: link.value || link.weight || 95
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

    const width = 1400;
    const height = 800;

    // Clear previous visualization
    d3.select(svgRef.current).selectAll('*').remove();

    const svg = d3.select(svgRef.current)
      .attr('width', width)
      .attr('height', height);

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
      return nodeIds.has(sourceId) && nodeIds.has(targetId) && link.value >= similarityScore;
    });

    // Create force simulation
    const simulation = d3.forceSimulation(filteredNodes)
      .force('link', d3.forceLink(filteredLinks).id((d: any) => d.id).distance(50))
      .force('charge', d3.forceManyBody().strength(-100))
      .force('center', d3.forceCenter(width / 2, height / 2))
      .force('collision', d3.forceCollide().radius(15));

    // Create links
    const link = g.append('g')
      .attr('class', 'links')
      .selectAll('line')
      .data(filteredLinks)
      .enter().append('line')
      .attr('stroke', '#999')
      .attr('stroke-opacity', 0.6)
      .attr('stroke-width', (d: any) => Math.sqrt(d.value - 70));

    // Create nodes
    const node = g.append('g')
      .attr('class', 'nodes')
      .selectAll('circle')
      .data(filteredNodes)
      .enter().append('circle')
      .attr('r', 8)
      .attr('fill', (d: Node) => {
        const colors: { [key: string]: string } = {
          'K': '#3b82f6',
          'Na': '#ef4444',
          'Ca': '#10b981',
          'IH': '#f59e0b',
          'KCa': '#8b5cf6',
          'Other': '#6b7280'
        };
        return colors[d.ion_class || 'Other'] || '#6b7280';
      })
      .attr('stroke', (d: Node) => d.supermodel ? '#fbbf24' : '#fff')
      .attr('stroke-width', (d: Node) => d.supermodel ? 3 : 1.5)
      .on('click', (event: any, d: any) => {
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

    // Simulation tick
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
  }, [networkData, ionClass, similarityScore, showICG, supermodel1, supermodel2, copiesNumber]);

  if (isLoading) {
    return (
      <div className="flex items-center justify-center min-h-screen bg-slate-50 dark:bg-slate-900">
        <div className="text-xl">Loading visualization...</div>
      </div>
    );
  }

  return (
    <div className="min-h-screen bg-slate-50 dark:bg-slate-900">
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

      <div className="flex">
        {/* Control Panel */}
        <div className="w-80 bg-white dark:bg-slate-800 border-r border-slate-200 dark:border-slate-700 h-[calc(100vh-57px)] overflow-y-auto">
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
                  <div className="w-4 h-4 bg-blue-500 rounded-full"></div>
                  <span className="text-slate-600 dark:text-slate-400">K Channel</span>
                </div>
                <div className="flex items-center space-x-2">
                  <div className="w-4 h-4 bg-red-500 rounded-full"></div>
                  <span className="text-slate-600 dark:text-slate-400">Na Channel</span>
                </div>
                <div className="flex items-center space-x-2">
                  <div className="w-4 h-4 bg-green-500 rounded-full"></div>
                  <span className="text-slate-600 dark:text-slate-400">Ca Channel</span>
                </div>
                <div className="flex items-center space-x-2">
                  <div className="w-4 h-4 bg-amber-500 rounded-full"></div>
                  <span className="text-slate-600 dark:text-slate-400">IH Channel</span>
                </div>
                <div className="flex items-center space-x-2">
                  <div className="w-4 h-4 bg-purple-500 rounded-full"></div>
                  <span className="text-slate-600 dark:text-slate-400">KCa Channel</span>
                </div>
                <div className="flex items-center space-x-2">
                  <div className="w-4 h-4 bg-gray-500 rounded-full"></div>
                  <span className="text-slate-600 dark:text-slate-400">Other</span>
                </div>
                <div className="flex items-center space-x-2">
                  <div className="w-4 h-4 border-2 border-amber-400 rounded-full"></div>
                  <span className="text-slate-600 dark:text-slate-400">Supermodel</span>
                </div>
              </div>
            </div>
          </div>
        </div>

        {/* Visualization Area */}
        <div className="flex-1 p-4">
          <div className="bg-white dark:bg-slate-800 rounded-lg shadow-lg p-4">
            <svg ref={svgRef}></svg>
          </div>
        </div>
      </div>
    </div>
  );
}
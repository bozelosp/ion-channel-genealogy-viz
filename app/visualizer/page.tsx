/* eslint-disable @typescript-eslint/no-explicit-any */
'use client';

import NextLink from 'next/link';
import { useEffect, useRef, useState, useMemo, useCallback, useDeferredValue } from 'react';
import * as d3 from 'd3';
import {
  formatIonClassForDisplay,
  getVisibleIonClassOptions,
  shouldSuppressIonClass,
  OTHER_ION_CLASS_LABEL,
} from '@/lib/ionClasses';

const ION_CLASS_NAMES = getVisibleIonClassOptions();
const ION_CLASS_FILTER_OPTIONS = ['all', ...ION_CLASS_NAMES];

interface Node {
  id: string;
  name?: string;
  label?: string;
  ion_class?: string;
  omnimodel1?: boolean;
  omnimodel2?: boolean;
  icg?: boolean;
  num_of_identicals?: number;
  original_model?: {
    ion_class?: string;
    ICG?: boolean;
    [key: string]: any;
  };
  identical_models?: Array<Record<string, unknown>>;
  omnimodel_location?: [channel: string, slug: string] | null;
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

type FixedLocationCircles = Record<number, Array<Array<[number, number, number]>>>;

const INITIAL_NETWORK_DATA_URL = '/data/network_data_initial.json';
const FULL_NETWORK_DATA_URL = '/data/network_data.json';
const INITIAL_MINIMUM_COPIES = 2;

function recordValue(value: unknown): Record<string, unknown> | null {
  return value !== null && typeof value === 'object' && !Array.isArray(value)
    ? value as Record<string, unknown>
    : null;
}

function datasetIdentifier(value: unknown): string | null {
  if (typeof value === 'string' && value.length > 0) return `string:${value}`;
  if (typeof value === 'number' && Number.isFinite(value)) return `number:${value}`;
  return null;
}

function validateNetworkDataset(value: unknown): { nodes: any[]; links: any[] } {
  const root = recordValue(value);
  const nodes = root?.nodes;
  const links = root?.links;
  if (!Array.isArray(nodes) || nodes.length === 0 || !Array.isArray(links) || links.length === 0) {
    throw new Error('Network dataset must contain non-empty node and link arrays');
  }

  const nodeIds = new Set<string>();
  const canonicalNodeIds = new Set<string>();
  const scientificIds = new Set<string>();
  for (const valueNode of nodes) {
    const node = recordValue(valueNode);
    const nodeId = datasetIdentifier(node?.id);
    const originalModel = recordValue(node?.original_model);
    const scientificId = originalModel?.unique_modelDB_mod_id;
    if (!node || !nodeId || typeof scientificId !== 'string' || scientificId.length === 0) {
      throw new Error('Network dataset contains an invalid node');
    }
    const canonicalNodeId = String(node.id);
    if (
      nodeIds.has(nodeId)
      || canonicalNodeIds.has(canonicalNodeId)
      || scientificIds.has(scientificId)
    ) {
      throw new Error('Network dataset contains duplicate node identifiers');
    }
    nodeIds.add(nodeId);
    canonicalNodeIds.add(canonicalNodeId);
    scientificIds.add(scientificId);
  }

  for (const valueLink of links) {
    const link = recordValue(valueLink);
    const source = datasetIdentifier(link?.source);
    const target = datasetIdentifier(link?.target);
    if (!link || !source || !target || !nodeIds.has(source) || !nodeIds.has(target)) {
      throw new Error('Network dataset contains a link with an invalid endpoint');
    }
  }

  return { nodes, links };
}

function processNetworkDataset(value: unknown): NetworkData {
  const rawData = validateNetworkDataset(value);
  return {
    nodes: rawData.nodes.map((node: any) => ({
      ...node,
      ion_class: node.ion_class || node.original_model?.ion_class || node.ion_channel_class || OTHER_ION_CLASS_LABEL,
      icg: typeof node.icg === 'boolean' ? node.icg : !!node.original_model?.ICG,
      omnimodel1: toBooleanFlag(node.omnimodel1 ?? node[LEGACY_MODEL_FIELD_1] ?? node.original_model?.[LEGACY_MODEL_KEY_1]),
      omnimodel2: toBooleanFlag(node.omnimodel2 ?? node[LEGACY_MODEL_FIELD_2] ?? node.original_model?.[LEGACY_MODEL_KEY_2]),
      num_of_identicals: node.num_of_identicals || 1,
    })),
    links: rawData.links.map((link: any) => ({
      ...link,
      weight: link.weight || link.value || 95,
    })),
  };
}

interface SelectionState {
  selectedIds: string[];
  sourceIds: string[];
  anchorId: string | null;
}

interface SourceViewState {
  node: Node;
  status: 'loading' | 'ready' | 'error';
  sourceCode?: string;
  error?: string;
}

function uniqueIdList(ids: Iterable<string>): string[] {
  const seen = new Set<string>();
  const result: string[] = [];
  for (const raw of ids) {
    const id = String(raw);
    if (!seen.has(id)) {
      seen.add(id);
      result.push(id);
    }
  }
  return result;
}

function normalizeSelectionState(state: SelectionState): SelectionState {
  const sourceIds = uniqueIdList(state.sourceIds);
  const sourceSet = new Set(sourceIds);
  const selectedIds = uniqueIdList(state.selectedIds.filter((id) => !sourceSet.has(id)));
  const anchorId = state.anchorId && (sourceSet.has(state.anchorId) || selectedIds.includes(state.anchorId))
    ? state.anchorId
    : (selectedIds[selectedIds.length - 1] ?? sourceIds[sourceIds.length - 1] ?? null);
  return { selectedIds, sourceIds, anchorId };
}

function selectionStatesEqual(a: SelectionState, b: SelectionState): boolean {
  if (a === b) return true;
  if (a.selectedIds.length !== b.selectedIds.length) return false;
  if (a.sourceIds.length !== b.sourceIds.length) return false;
  if (a.anchorId !== b.anchorId) return false;
  for (let i = 0; i < a.selectedIds.length; i += 1) {
    if (a.selectedIds[i] !== b.selectedIds[i]) return false;
  }
  for (let i = 0; i < a.sourceIds.length; i += 1) {
    if (a.sourceIds[i] !== b.sourceIds[i]) return false;
  }
  return true;
}

// Calculate node radius based on num_of_identicals (from original implementation)
function calculateNodeRadius(node: Node): number {
  const numIdenticals = node.num_of_identicals || 1;
  const radius = (0.2375 * Math.log(numIdenticals) / Math.log(1.09) + 1.325 + 
                  0.3925 * Math.log(numIdenticals) / Math.log(1.35) + 3) / 2;
  return radius;
}

// Group splitting constants (matching original)
const FORCE_STRENGTH = 0.1675;

const LEGACY_MODEL_PREFIX = 'Super' + 'model';
const LEGACY_MODEL_KEY_1 = `${LEGACY_MODEL_PREFIX} 1`;
const LEGACY_MODEL_KEY_2 = `${LEGACY_MODEL_PREFIX} 2`;
const LEGACY_MODEL_FIELD_1 = (LEGACY_MODEL_PREFIX + '1').toLowerCase();
const LEGACY_MODEL_FIELD_2 = (LEGACY_MODEL_PREFIX + '2').toLowerCase();

const toBooleanFlag = (value: unknown): boolean => (typeof value === 'boolean' ? value : !!value);

const extractModelFlag = (
  node: Node,
  key: 'omnimodel1' | 'omnimodel2',
  legacyField: string,
  legacyKey: string
): boolean => {
  const rawNode = node as unknown as Record<string, unknown>;
  const primaryValue = rawNode[key];
  const fallbackValue = rawNode[legacyField];
  const legacyValue = node.original_model?.[legacyKey];
  return toBooleanFlag(primaryValue ?? fallbackValue ?? legacyValue);
};

const OMNIMODEL_BASE_CHANNELS = new Set(['Na', 'K', 'Ca', 'IH', 'KCa']);

const OMNIMODEL_CHANNEL_ALIASES: Record<string, string> = {
  na: 'Na',
  k: 'K',
  ca: 'Ca',
  ih: 'IH',
  h: 'IH',
  kca: 'KCa',
};

function normaliseChannel(raw?: string | null): string | null {
  if (!raw) return null;
  const trimmed = raw.trim();
  if (!trimmed) return null;
  const alias = OMNIMODEL_CHANNEL_ALIASES[trimmed.toLowerCase()];
  const resolved = alias ?? trimmed;
  if (!OMNIMODEL_BASE_CHANNELS.has(resolved)) {
    return null;
  }
  return resolved;
}

interface OmnimodelLocation {
  channel: string;
  slug: string;
}

function extractLocationFromIcgPath(pathValue: unknown): OmnimodelLocation | null {
  if (typeof pathValue !== 'string') return null;
  const segments = pathValue.split('/').filter(Boolean);
  if (segments.length < 2) return null;
  const channelSegment = segments.find((segment) => segment.startsWith('icg-channels-'));
  if (!channelSegment) return null;
  const channel = normaliseChannel(channelSegment.replace('icg-channels-', ''));
  if (!channel) return null;
  const channelIndex = segments.indexOf(channelSegment);
  const slugSegment = segments[channelIndex + 1];
  if (!slugSegment) return null;
  return { channel, slug: slugSegment };
}

function deriveFallbackLocation(node: Node | null | undefined): OmnimodelLocation | null {
  if (!node?.original_model) return null;
  const channel = normaliseChannel(node.original_model.ion_class ?? node.ion_class);
  if (!channel) return null;
  const modelDbDir = node.original_model.modelDB_dir as string | undefined;
  const modFilename = node.original_model.mod_filename as string | undefined;
  if (!modelDbDir || !modFilename) return null;
  const filename = modFilename.endsWith('.mod') ? modFilename : `${modFilename}.mod`;
  const slugBase = `${modelDbDir}_${filename}`;
  return { channel, slug: slugBase };
}

function deriveOmnimodelLocation(node: Node | null | undefined): OmnimodelLocation | null {
  if (!node) return null;
  const compactLocation = node.omnimodel_location;
  if (Array.isArray(compactLocation) && compactLocation.length === 2) {
    const channel = normaliseChannel(compactLocation[0]);
    const slug = compactLocation[1];
    if (channel && typeof slug === 'string' && slug.length > 0) return { channel, slug };
  }
  const icgEntries = (node.original_model as any)?.ICG_entries;
  if (Array.isArray(icgEntries)) {
    for (const entry of icgEntries) {
      const info = entry?.info;
      const loc = extractLocationFromIcgPath(info?.mod_filepath ?? info?.modFilename ?? null);
      if (loc) {
        return loc;
      }
      const locFromFilename = extractLocationFromIcgPath(entry?.mod_filepath ?? entry?.path);
      if (locFromFilename) {
        return locFromFilename;
      }
    }
  }
  return deriveFallbackLocation(node);
}

function makeOmnimodelKey(location: OmnimodelLocation): string {
  return `${location.channel}::${location.slug}`;
}

function getOmnimodelPath(
  node: Node | null | undefined,
  availableKeys?: Set<string>,
): string | null {
  const location = deriveOmnimodelLocation(node);
  if (!location) return null;
  if (availableKeys && !availableKeys.has(makeOmnimodelKey(location))) {
    return null;
  }
  return `/omnimodels/${encodeURIComponent(location.channel)}/${encodeURIComponent(location.slug)}`;
}

function areSetsEqual<T>(a: Set<T> | null | undefined, b: Set<T>): boolean {
  if (!a || a.size !== b.size) {
    return false;
  }
  for (const item of b) {
    if (!a.has(item)) {
      return false;
    }
  }
  return true;
}

function getNodeUniqueId(node: Node | null | undefined): string | undefined {
  return node?.original_model?.unique_modelDB_mod_id as string | undefined;
}

function getNodeDisplayName(node: Node | null | undefined): string {
  if (!node) return '';
  if (typeof node.id === 'string' && node.id.trim().length > 0) {
    return node.id;
  }
  const uniqueId = getNodeUniqueId(node);
  if (uniqueId) {
    return uniqueId;
  }
  return 'Unknown node';
}

function parseNumeric(value: unknown): number | null {
  if (typeof value === 'number' && Number.isFinite(value)) return value;
  if (typeof value === 'string') {
    const parsed = parseInt(value, 10);
    if (Number.isFinite(parsed)) return parsed;
  }
  return null;
}

// Normalize free-text queries: trim and right-strip a trailing .mod
function normalizeQueryTerm(raw: string): string {
  const trimmed = (raw || '').trim();
  let lower = trimmed.toLowerCase();
  if (lower.endsWith('.mod')) {
    lower = lower.slice(0, -4);
  }
  return lower;
}

function getCandidateYear(node: Node, fallback?: Record<string, unknown>): number {
  const yearFromNode = parseNumeric(node.original_model?.Year);
  if (yearFromNode !== null) return yearFromNode;
  const yearFromFallback = parseNumeric(fallback?.Year);
  if (yearFromFallback !== null) return yearFromFallback;
  return Number.POSITIVE_INFINITY;
}

function getCandidateModelDbDir(node: Node, fallback?: Record<string, unknown>): number {
  const dirFromNode = parseNumeric(node.original_model?.modelDB_dir);
  if (dirFromNode !== null) return dirFromNode;
  const dirFromFallback = parseNumeric(fallback?.modelDB_dir);
  if (dirFromFallback !== null) return dirFromFallback;
  return Number.POSITIVE_INFINITY;
}

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

// Get unique filter combinations based on active filters (matching original)
function getUniqueFilterCombinations(selectedIonClasses: Set<string>, showICG: boolean, omnimodel1: boolean, omnimodel2: boolean): string[][] {
  // Get active filters (only if actually active)
  const activeFilters: string[] = [];
  if (omnimodel1) activeFilters.push('Omnimodel 1');
  if (omnimodel2) activeFilters.push('Omnimodel 2');
  if (showICG) activeFilters.push('ICG entry');

  // Get active ion classes
  let activeIonClasses = selectedIonClasses.has('all')
    ? ['All']
    : Array.from(selectedIonClasses).filter((cls) => !shouldSuppressIonClass(cls));

  if (!selectedIonClasses.has('all') && activeIonClasses.length === 0) {
    activeIonClasses = ['All'];
  }

  // Initialize list for unique combinations
  let uniqueFilterCombinations: string[][] = [];

  // Generate all unique combinations of active filters
  const filterCombinations = getCombinations(activeFilters);

  // Combine each filter combination with each ion class (matching original)
  for (const ionClass of activeIonClasses) {
    for (const filterCombo of filterCombinations) {
      uniqueFilterCombinations.push([ionClass, ...filterCombo]);
    }
  }

  // Add single ion classes as their own group (this is the key from original line 43!)
  uniqueFilterCombinations = [...uniqueFilterCombinations, ...activeIonClasses.map(c => [c])];

  // Remove duplicates if any
  const seen = new Set();
  uniqueFilterCombinations = uniqueFilterCombinations.filter(combo => {
    const key = combo.join(',');
    if (seen.has(key)) return false;
    seen.add(key);
    return true;
  });

  return uniqueFilterCombinations;
}

// Assign nodes to groups based on their properties
function assignNodesToGroups(
  nodes: Node[], 
  uniqueFilterCombinations: string[][], 
  selectedIonClasses: Set<string>,
  showICG: boolean,
  omnimodel1: boolean,
  omnimodel2: boolean
): { [key: string]: Node[] } {
  const nodesGroupedByFilter: { [key: string]: Node[] } = {};

  uniqueFilterCombinations.forEach(filterCombination => {
    const groupKey = filterCombination.join(',');

    nodes.forEach(node => {
      // Use original_model fields consistently (matching original implementation)
      const nodeIonClass = node.original_model?.ion_class || node.ion_class;
      if (shouldSuppressIonClass(nodeIonClass)) {
        return;
      }
      const nodeOm1 = extractModelFlag(node, 'omnimodel1', LEGACY_MODEL_FIELD_1, LEGACY_MODEL_KEY_1);
      const nodeOm2 = extractModelFlag(node, 'omnimodel2', LEGACY_MODEL_FIELD_2, LEGACY_MODEL_KEY_2);
      const nodeIcg = !!(node.icg || node.original_model?.ICG);

      // Match original logic: filter_state['All']?.filter_value || filter_combination.includes(ion_class)
      const ionClassMatch = selectedIonClasses.has('all') || filterCombination.includes(nodeIonClass || '');

      let status = true;

      if (ionClassMatch) {
        // Check Omnimodel 1
        if (omnimodel1) {
          if (filterCombination.includes('Omnimodel 1')) {
            if (!nodeOm1) status = false;
          } else {
            if (nodeOm1) status = false;
          }
        }

        // Check Omnimodel 2
        if (omnimodel2) {
          if (filterCombination.includes('Omnimodel 2')) {
            if (!nodeOm2) status = false;
          } else {
            if (nodeOm2) status = false;
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

function filterNetworkNodes(
  nodes: Node[],
  selectedIonClasses: Set<string>,
  showICG: boolean,
  minimumCopies: number,
): Node[] {
  return nodes.filter((node) => {
    const nodeIonClass = node.original_model?.ion_class || node.ion_class;
    if (shouldSuppressIonClass(nodeIonClass)) return false;
    if (!selectedIonClasses.has('all') && (!nodeIonClass || !selectedIonClasses.has(nodeIonClass))) {
      return false;
    }
    if (showICG && !(node.icg || node.original_model?.ICG)) return false;
    return (node.num_of_identicals ?? 0) >= minimumCopies;
  });
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
  let myCirclePositions: any = null;
  
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
    return groups.map(() => [0.5, 0.5]);
  }

  return myCirclePositions.map((d: any) => [avoidBorders(d[1]), avoidBorders(d[2])]);
}

// Sort and position groups (matching original implementation)
function sortAndPositionGroups(
  nodesGroupedByFilter: { [key: string]: Node[] },
  fixedLocationCircles: any
): { [key: string]: [number, number] } {
  // If only one group, return empty object (no positioning needed)
  const numberOfGroups = Object.keys(nodesGroupedByFilter).length;
  if (numberOfGroups <= 1) {
    return {};
  }

  const groups: any[] = [];
  
  for (const key in nodesGroupedByFilter) {
    if (nodesGroupedByFilter[key].length > 0) {
      const groupRadiusValues = nodesGroupedByFilter[key].map(node => calculateNodeRadius(node));
      const totalGroupRadius = groupRadiusValues.reduce((a, b) => a + b, 0);
      groups.push([key, Math.pow(totalGroupRadius, 0.8), nodesGroupedByFilter[key]]);
    }
  }

  // Sort groups by size (ascending)
  groups.sort((a, b) => a[1] - b[1]);

  // Position the groups using fixed locations
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
// Limit concurrency for network requests
async function runLimited<T>(items: T[], limit: number, worker: (item: T) => Promise<void>): Promise<void> {
  let index = 0;
  const workers = Array.from({ length: Math.min(limit, items.length) }, async () => {
    while (true) {
      const i = index++;
      if (i >= items.length) break;
      await worker(items[i]);
    }
  });
  await Promise.all(workers);
}

async function fetchSourceCode(
  nodeIds: string[],
  networkData: NetworkData,
  fetchedFiles: { [key: string]: { source_code: string } }
): Promise<{
  files: { [key: string]: { source_code: string } },
  failures: { nodeId: string, uniqueId?: string, name?: string, error?: string }[]
}> {
  const updatedFetchedFiles = { ...fetchedFiles };
  const failures: { nodeId: string, uniqueId?: string, name?: string, error?: string }[] = [];
  const nodeById = new Map(networkData.nodes.map((node) => [String(node.id), node]));

  const uniqueNodeIds = Array.from(new Set(nodeIds));
  const concurrency = 8;

  await runLimited(uniqueNodeIds, concurrency, async (rawId) => {
    try {
      const nodeId = String(rawId);
      if (updatedFetchedFiles[nodeId]) return;
      const nodeData = nodeById.get(nodeId);
      if (!nodeData) {
        console.warn(`Node with ID ${nodeId} not found.`);
        failures.push({ nodeId, error: 'Node not found' });
        return;
      }
      const uniqueId = (nodeData as any).original_model?.unique_modelDB_mod_id;
      if (!uniqueId) {
        console.warn(`No unique_modelDB_mod_id found for node ${nodeId}`);
        failures.push({ nodeId, name: (nodeData as any).name || nodeData.id, error: 'No unique_modelDB_mod_id' });
        return;
      }
      // Use path param form so the CDN forwards it to the origin
      const response = await fetch(`/api/source-code/${encodeURIComponent(uniqueId)}`);
      if (!response.ok) {
        throw new Error(`HTTP ${response.status}`);
      }
      const sourceCode = await response.text();
      updatedFetchedFiles[nodeId] = { source_code: sourceCode };
    } catch (error) {
      const nodeId = String(rawId);
      console.error(`Error fetching file for node ${nodeId}:`, error);
      const nodeData = nodeById.get(nodeId);
      failures.push({
        nodeId,
        name: (nodeData as any)?.name || nodeId,
        uniqueId: (nodeData as any)?.original_model?.unique_modelDB_mod_id,
        error: error instanceof Error ? error.message : 'Unknown error'
      });
    }
  });

  return { files: updatedFetchedFiles, failures };
}

// Generate bounded diff HTML away from the main UI thread.
type DiffWorkerReply = { id: number; html?: string; error?: string };
type PendingDiff = {
  resolve: (html: string) => void;
  timeout: ReturnType<typeof setTimeout>;
};

const DIFF_WORKER_TIMEOUT_MS = 8_000;
const MAX_DIFF_SOURCE_CHARACTERS = 64 * 1024;
let diffWorker: Worker | null = null;
let nextDiffRequestId = 1;
const pendingDiffs = new Map<number, PendingDiff>();

function diffErrorDocument(message: string): string {
  const escaped = message.replace(/[&<>"']/g, (character) => ({
    '&': '&amp;',
    '<': '&lt;',
    '>': '&gt;',
    '"': '&quot;',
    "'": '&#039;',
  })[character] as string);
  return `<!doctype html><html><body><div style="padding:16px;border:1px solid #fecaca;border-radius:8px;background:#fef2f2;color:#b91c1c;font-family:system-ui,sans-serif">Diff unavailable: ${escaped}</div></body></html>`;
}

function diffPlaceholderDocument(message: string, theme: 'light' | 'dark'): string {
  const escaped = message.replace(/[&<>"']/g, (character) => ({
    '&': '&amp;',
    '<': '&lt;',
    '>': '&gt;',
    '"': '&quot;',
    "'": '&#039;',
  })[character] as string);
  return `<!doctype html><html data-theme="${theme}"><head><meta charset="utf-8"/><meta name="viewport" content="width=device-width,initial-scale=1"/><style>:root{--bg:#fff;--text:#0f172a;--border:#e2e8f0}:root[data-theme='dark']{--bg:#0f172a;--text:#e2e8f0;--border:#334155}html,body{margin:0;padding:16px;background:var(--bg);color:var(--text);font-family:system-ui,sans-serif}.card{border:1px solid var(--border);border-radius:8px;padding:16px}.title{font-weight:700;margin-bottom:8px}.muted{opacity:.8}</style></head><body><div class="card"><div class="title">Diff unavailable</div><div class="muted">${escaped}</div></div></body></html>`;
}

function resetDiffWorker(message: string): void {
  const errorHtml = diffErrorDocument(message);
  for (const { resolve, timeout } of pendingDiffs.values()) {
    clearTimeout(timeout);
    resolve(errorHtml);
  }
  pendingDiffs.clear();
  diffWorker?.terminate();
  diffWorker = null;
}

function getDiffWorker(): Worker {
  if (diffWorker) return diffWorker;

  const worker = new Worker(new URL('./diff.worker.ts', import.meta.url), { type: 'module' });
  worker.onmessage = (event: MessageEvent<DiffWorkerReply>) => {
    const pending = pendingDiffs.get(event.data.id);
    if (!pending) return;

    clearTimeout(pending.timeout);
    pendingDiffs.delete(event.data.id);
    pending.resolve(
      event.data.html ?? diffErrorDocument(event.data.error ?? 'Diff generation failed'),
    );
  };
  worker.onerror = () => resetDiffWorker('Diff worker failed');
  diffWorker = worker;
  return worker;
}

async function generateDiffHtml(
  sourceCode: string,
  targetCode: string,
  theme: 'light' | 'dark' = 'light',
): Promise<string> {
  if (
    sourceCode.length > MAX_DIFF_SOURCE_CHARACTERS ||
    targetCode.length > MAX_DIFF_SOURCE_CHARACTERS
  ) {
    return diffErrorDocument('A source file exceeds the 64 KiB comparison limit');
  }

  try {
    const worker = getDiffWorker();
    const id = nextDiffRequestId;
    nextDiffRequestId += 1;

    return await new Promise<string>((resolve) => {
      const timeout = setTimeout(() => {
        if (pendingDiffs.has(id)) {
          resetDiffWorker('Diff generation exceeded the 8-second limit');
        }
      }, DIFF_WORKER_TIMEOUT_MS);

      pendingDiffs.set(id, { resolve, timeout });
      worker.postMessage({ id, source: sourceCode, target: targetCode, theme });
    });
  } catch (error) {
    return diffErrorDocument(error instanceof Error ? error.message : 'Diff generation failed');
  }
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

export default function Visualizer() {
  // Hard safety limits to keep UI responsive and servers happy
  const MAX_SELECTED_NODES = 50;
  const MAX_COMBINATIONS = 200;
  const svgRef = useRef<SVGSVGElement>(null);
  const infoFrozenRef = useRef<boolean>(false);
  const sourceDialogRef = useRef<HTMLDivElement>(null);
  const sourceDialogCloseRef = useRef<HTMLButtonElement>(null);
  const sourceDialogPreviousFocusRef = useRef<HTMLElement | null>(null);
  const sourceFetchRequestIdRef = useRef(0);
  const [networkData, setNetworkData] = useState<NetworkData | null>(null);
  const [hasFullNetworkData, setHasFullNetworkData] = useState(false);
  const [isFullNetworkLoading, setIsFullNetworkLoading] = useState(false);
  const [fullNetworkError, setFullNetworkError] = useState<string | null>(null);
  const [fixedLocationCircles, setFixedLocationCircles] = useState<FixedLocationCircles>({});
  const [isGroupedLayoutLoading, setIsGroupedLayoutLoading] = useState(false);
  const [groupedLayoutError, setGroupedLayoutError] = useState<string | null>(null);
  const [layoutRetry, setLayoutRetry] = useState(0);
  const [selectedIonClasses, setSelectedIonClasses] = useState<Set<string>>(new Set(['all']));
  const [similarityScore, setSimilarityScore] = useState<number>(95);
  const deferredSimilarityScore = useDeferredValue(similarityScore);
  // Persist similarity in ?sim= and read initial value from URL
  useEffect(() => {
    try {
      const url = new URL(window.location.href);
      const simParam = url.searchParams.get('sim');
      if (simParam) {
        const num = parseInt(simParam, 10);
        if (!Number.isNaN(num) && num >= 75 && num <= 100) {
          setSimilarityScore(num);
        }
      }
    } catch {}
  }, []);
  useEffect(() => {
    try {
      const url = new URL(window.location.href);
      url.searchParams.set('sim', String(similarityScore));
      window.history.replaceState({}, '', url.toString());
    } catch {}
  }, [similarityScore]);
  const [copiesNumber, setCopiesNumber] = useState<number>(2);
  const effectiveCopiesNumber = !hasFullNetworkData && copiesNumber < INITIAL_MINIMUM_COPIES
    ? INITIAL_MINIMUM_COPIES
    : copiesNumber;
  const deferredCopiesNumber = useDeferredValue(effectiveCopiesNumber);
  const [showICG, setShowICG] = useState<boolean>(false);
  const [omnimodel1, setOmnimodel1] = useState<boolean>(false);
  const [omnimodel2, setOmnimodel2] = useState<boolean>(false);
  const [selectedNode, setSelectedNode] = useState<Node | null>(null);
  const [isLoading, setIsLoading] = useState<boolean>(true);
  const [dataLoadError, setDataLoadError] = useState<string | null>(null);
  const [isMobile, setIsMobile] = useState<boolean>(false);
  const [containerDimensions, setContainerDimensions] = useState<{width: number, height: number}>({width: 0, height: 0});
  const [searchTerm, setSearchTerm] = useState<string>('');
  const [searchMatches, setSearchMatches] = useState<Set<string>>(new Set());
  const [visibleNodeIds, setVisibleNodeIds] = useState<Set<string>>(new Set());
  const [selectionState, setSelectionState] = useState<SelectionState>({ selectedIds: [], sourceIds: [], anchorId: null });
  const selectionStateRef = useRef<SelectionState>(selectionState);
  
  // Detect system dark mode preference
  const [resolvedTheme, setResolvedTheme] = useState<'light' | 'dark'>('light');
  useEffect(() => {
    const mq = window.matchMedia('(prefers-color-scheme: dark)');
    setResolvedTheme(mq.matches ? 'dark' : 'light');
    const handler = (e: MediaQueryListEvent) => setResolvedTheme(e.matches ? 'dark' : 'light');
    mq.addEventListener('change', handler);
    return () => mq.removeEventListener('change', handler);
  }, []);
  
  // Source code comparison states
  const fetchedFilesRef = useRef<{[key: string]: {source_code: string}}>({});
  const [showDiffView, setShowDiffView] = useState<boolean>(false);
  const [sourceView, setSourceView] = useState<SourceViewState | null>(null);
  // Failed source fetches banner state
  const [failedNodes, setFailedNodes] = useState<{nodeId: string, uniqueId?: string, name?: string, error?: string}[]>([]);
  const [showFailureBanner, setShowFailureBanner] = useState<boolean>(true);
  const [isFailureExpanded, setIsFailureExpanded] = useState<boolean>(false);
  const [retryingFailed, setRetryingFailed] = useState<boolean>(false);

  // View control
  const [fitToViewFunction, setFitToViewFunction] = useState<(() => void) | null>(null);
  
  // Context menu state
  const [contextMenu, setContextMenu] = useState<{x: number, y: number, node: Node} | null>(null);
  const [diffCombinations, setDiffCombinations] = useState<{source: Node, target: Node, html?: string}[]>([]);
  const [currentCombinationIndex, setCurrentCombinationIndex] = useState<number>(0);
  const [isGeneratingDiffs, setIsGeneratingDiffs] = useState<boolean>(false);
  const inFlightDiffKeysRef = useRef(new Set<string>());
  const [isInfoFrozen, setIsInfoFrozen] = useState<boolean>(false);
  const [leftFontScale, setLeftFontScale] = useState<number>(1.05);
  // Subgraph selection function (set by D3 setup, used by keyboard 'A')
  const [subgraphSelectFn, setSubgraphSelectFn] = useState<((n: Node) => Promise<void> | void) | null>(null);
  // Selection summary (right panel)
  const [selectionSummary, setSelectionSummary] = useState<{
    nodes: number;
    copies: number;
    rows: Array<{ node_id: string; unique_id: string; name?: string; ion_class?: string; copies: number }>;
  }>({ nodes: 0, copies: 0, rows: [] });
  // Ref to clear current marquee rectangle (for Escape)
  const clearBrushRef = useRef<(() => void) | null>(null);
  const isBrushingRef = useRef<boolean>(false);
  const modifierStateRef = useRef({ meta: false, ctrl: false, shift: false, alt: false });
  const marqueeStateRef = useRef<{ active: boolean; startX: number; startY: number; mode: 'replace' | 'add' | 'subtract' }>({ active: false, startX: 0, startY: 0, mode: 'replace' });
  const marqueeRectRef = useRef<d3.Selection<SVGRectElement, unknown, null, undefined> | null>(null);
  
  // Group summary state
  const [groupSummaries, setGroupSummaries] = useState<{key: string, nodeCount: number}[]>([]);

  const storeFetchedFiles = useCallback((files: {[key: string]: {source_code: string}}) => {
    fetchedFilesRef.current = files;
  }, []);

  const closeSourceView = useCallback(() => {
    sourceFetchRequestIdRef.current += 1;
    setSourceView(null);
    const previousFocus = sourceDialogPreviousFocusRef.current;
    sourceDialogPreviousFocusRef.current = null;
    window.requestAnimationFrame(() => previousFocus?.focus());
  }, []);

  const openSourceForNode = useCallback(async (node: Node | null | undefined) => {
    if (!node || !networkData) return;

    const nodeId = String(node.id);
    const requestId = sourceFetchRequestIdRef.current + 1;
    sourceFetchRequestIdRef.current = requestId;
    const activeElement = document.activeElement instanceof HTMLElement ? document.activeElement : null;
    if (!activeElement || !sourceDialogRef.current?.contains(activeElement)) {
      sourceDialogPreviousFocusRef.current = activeElement;
    }
    setContextMenu(null);
    setSourceView({ node, status: 'loading' });

    try {
      const { files, failures } = await fetchSourceCode([nodeId], networkData, fetchedFilesRef.current);
      if (sourceFetchRequestIdRef.current !== requestId) return;

      storeFetchedFiles(files);
      if (failures.length > 0) {
        setFailedNodes((previous) => [
          ...previous.filter((failure) => failure.nodeId !== nodeId),
          ...failures,
        ]);
        setShowFailureBanner(true);
        setSourceView({
          node,
          status: 'error',
          error: failures[0]?.error || 'Source code could not be loaded.',
        });
        return;
      }

      setFailedNodes((previous) => previous.filter((failure) => failure.nodeId !== nodeId));
      const sourceCode = files[nodeId]?.source_code;
      if (typeof sourceCode !== 'string') {
        setSourceView({ node, status: 'error', error: 'Source code was not available for this node.' });
        return;
      }

      setSourceView({ node, status: 'ready', sourceCode });
    } catch (error) {
      if (sourceFetchRequestIdRef.current !== requestId) return;
      setSourceView({
        node,
        status: 'error',
        error: error instanceof Error ? error.message : 'Source code could not be loaded.',
      });
    }
  }, [networkData, storeFetchedFiles]);

  const sourceViewNodeId = sourceView ? String(sourceView.node.id) : null;
  useEffect(() => {
    if (!sourceViewNodeId) return;
    const frame = window.requestAnimationFrame(() => sourceDialogCloseRef.current?.focus());
    return () => window.cancelAnimationFrame(frame);
  }, [sourceViewNodeId]);

  const needsGroupedLayout = omnimodel1 || omnimodel2 || showICG || selectedIonClasses.size > 1;
  const groupedLayoutCount = useMemo(() => {
    if (!networkData || !needsGroupedLayout) return 0;
    const nodes = filterNetworkNodes(
      networkData.nodes,
      selectedIonClasses,
      showICG,
      deferredCopiesNumber,
    );
    const combinations = getUniqueFilterCombinations(
      selectedIonClasses,
      showICG,
      omnimodel1,
      omnimodel2,
    );
    return Object.keys(assignNodesToGroups(
      nodes,
      combinations,
      selectedIonClasses,
      showICG,
      omnimodel1,
      omnimodel2,
    )).length;
  }, [
    deferredCopiesNumber,
    needsGroupedLayout,
    networkData,
    omnimodel1,
    omnimodel2,
    selectedIonClasses,
    showICG,
  ]);

  useEffect(() => {
    setGroupedLayoutError(null);
    if (groupedLayoutCount <= 1 || fixedLocationCircles[groupedLayoutCount]) {
      setIsGroupedLayoutLoading(false);
      return;
    }

    const controller = new AbortController();
    setIsGroupedLayoutLoading(true);
    fetch(`/data/fixed-location-circles/${groupedLayoutCount}.json`, { signal: controller.signal })
      .then(async (response) => {
        if (!response.ok) throw new Error(`HTTP ${response.status}`);
        const candidates = await response.json();
        if (
          !Array.isArray(candidates)
          || candidates.length === 0
          || candidates.some((candidate) => (
            !Array.isArray(candidate)
            || candidate.length !== groupedLayoutCount
            || candidate.some((circle) => (
              !Array.isArray(circle)
              || circle.length !== 3
              || circle.some((value) => typeof value !== 'number' || !Number.isFinite(value))
            ))
          ))
        ) {
          throw new Error(`Unexpected layout data for ${groupedLayoutCount} groups`);
        }
        if (controller.signal.aborted) return;
        setFixedLocationCircles((current) => ({
          ...current,
          [groupedLayoutCount]: candidates,
        }));
      })
      .catch((error) => {
        if (controller.signal.aborted) return;
        console.warn('Failed to load grouped layout positions:', error);
        setGroupedLayoutError('The grouped layout could not be loaded. The graph is hidden until its positions are available.');
      })
      .finally(() => {
        if (!controller.signal.aborted) setIsGroupedLayoutLoading(false);
      });
    return () => controller.abort();
  }, [fixedLocationCircles, groupedLayoutCount, layoutRetry]);
  
  useEffect(() => {
    const updateModifierState = (event: KeyboardEvent) => {
      modifierStateRef.current = {
        meta: event.metaKey,
        ctrl: event.ctrlKey,
        shift: event.shiftKey,
        alt: event.altKey,
      };
    };

    const handleKeyDown = (event: KeyboardEvent) => updateModifierState(event);
    const handleKeyUp = (event: KeyboardEvent) => updateModifierState(event);
    const handleBlur = () => {
      modifierStateRef.current = { meta: false, ctrl: false, shift: false, alt: false };
    };

    window.addEventListener('keydown', handleKeyDown, true);
    window.addEventListener('keyup', handleKeyUp, true);
    window.addEventListener('blur', handleBlur);
    return () => {
      window.removeEventListener('keydown', handleKeyDown, true);
      window.removeEventListener('keyup', handleKeyUp, true);
      window.removeEventListener('blur', handleBlur);
    };
  }, []);

  const selectedNodeOmnimodel1 = selectedNode ? extractModelFlag(selectedNode, 'omnimodel1', LEGACY_MODEL_FIELD_1, LEGACY_MODEL_KEY_1) : false;
  const selectedNodeOmnimodel2 = selectedNode ? extractModelFlag(selectedNode, 'omnimodel2', LEGACY_MODEL_FIELD_2, LEGACY_MODEL_KEY_2) : false;
  const selectedNodeAuthors = (() => {
    const raw = selectedNode?.original_model?.Authors;
    if (typeof raw !== 'string') return '';
    const trimmed = raw.trim();
    return trimmed.length > 0 ? trimmed : '';
  })();
  const [availableOmnimodelKeys, setAvailableOmnimodelKeys] = useState<Set<string> | null>(null);
  const selectedNodeDisplayName = useMemo(() => getNodeDisplayName(selectedNode), [selectedNode]);
  const currentDiffCombination = diffCombinations[currentCombinationIndex] ?? null;
  const currentSourceIonClass = currentDiffCombination?.source
    ? formatIonClassForDisplay(currentDiffCombination.source.original_model?.ion_class || currentDiffCombination.source.ion_class)
    : null;
  const currentTargetIonClass = currentDiffCombination?.target
    ? formatIonClassForDisplay(currentDiffCombination.target.original_model?.ion_class || currentDiffCombination.target.ion_class)
    : null;

  const syncSelectionToDom = useCallback((state: SelectionState) => {
    if (!svgRef.current) return;
    const svg = d3.select(svgRef.current);
    const baseFill = resolvedTheme === 'dark' ? '#3b82f6' : '#00BFFF';
    const baseStroke = resolvedTheme === 'dark' ? '#475569' : '#aaa';
    const selectedFill = resolvedTheme === 'dark' ? '#3b82f6' : '#00BFFF';
    const selectedStroke = resolvedTheme === 'dark' ? '#1e40af' : '#215885';
    const sourceFill = resolvedTheme === 'dark' ? '#fbbf24' : '#ffd700';
    const sourceStroke = resolvedTheme === 'dark' ? '#1e40af' : '#215885';

    const selectedSet = new Set(state.selectedIds);
    const sourceSet = new Set(state.sourceIds);

    svg.selectAll<SVGCircleElement, Node>('.graph-node').each(function (d) {
      const circle = d3.select<SVGCircleElement, Node>(this);
      const nodeId = String(d.id);
      const isSource = sourceSet.has(nodeId);
      const isSelected = selectedSet.has(nodeId);
      circle.classed('source-node', isSource);
      circle.classed('selected-node', isSelected && !isSource);

      let fill = baseFill;
      let stroke = baseStroke;
      let strokeWidth = 1;

      if (isSource) {
        fill = sourceFill;
        stroke = sourceStroke;
        strokeWidth = 2;
      } else if (isSelected) {
        fill = selectedFill;
        stroke = selectedStroke;
        strokeWidth = 2;
      }

      if (!circle.classed('search-override') || isSource) {
        circle.attr('fill', fill);
      }
      circle.attr('stroke', stroke).attr('stroke-width', strokeWidth);
    });
  }, [resolvedTheme]);

  const commitSelectionState = useCallback((updater: (prev: SelectionState) => SelectionState) => {
    setSelectionState((prev) => {
      const updated = normalizeSelectionState(updater(prev));
      if (selectionStatesEqual(prev, updated)) {
        return prev;
      }
      return updated;
    });
  }, []);

  useEffect(() => {
    selectionStateRef.current = selectionState;
    syncSelectionToDom(selectionState);
  }, [selectionState, syncSelectionToDom]);

  const nodeById = useMemo(() => {
    const map = new Map<string, Node>();
    if (networkData) {
      networkData.nodes.forEach((node) => {
        map.set(String(node.id), node);
      });
    }
    return map;
  }, [networkData]);

  const selectionBuckets = useMemo(() => {
    const sources: Node[] = [];
    const selected: Node[] = [];
    const all: Node[] = [];
    if (nodeById.size === 0) {
      return { sources, selected, all };
    }
    selectionState.sourceIds.forEach((id) => {
      const node = nodeById.get(id);
      if (node) {
        sources.push(node);
        all.push(node);
      }
    });
    selectionState.selectedIds.forEach((id) => {
      const node = nodeById.get(id);
      if (node && !sources.some((n) => n.id === node.id)) {
        selected.push(node);
        all.push(node);
      }
    });
    return { sources, selected, all };
  }, [nodeById, selectionState]);

  useEffect(() => {
    if (selectionBuckets.all.length === 0) {
      setSelectionSummary({ nodes: 0, copies: 0, rows: [] });
      return;
    }

    const rows: Array<{ node_id: string; unique_id: string; name?: string; ion_class?: string; copies: number }> = [];
    let copiesTotal = 0;
    selectionBuckets.all.forEach((node) => {
      const node_id = String(node.id);
      const unique_id = (node.original_model?.unique_modelDB_mod_id as string | undefined) || '';
      const name = (node.name as string | undefined) || (node.label as string | undefined);
      const ion_class = formatIonClassForDisplay(node.original_model?.ion_class || node.ion_class || '');
      const copies = (node.num_of_identicals ?? 1) || 1;
      copiesTotal += copies;
      rows.push({ node_id, unique_id, name, ion_class, copies });
    });

    setSelectionSummary({ nodes: rows.length, copies: copiesTotal, rows });
  }, [selectionBuckets]);

  const selectedNodeOmnimodelHref = availableOmnimodelKeys
    ? getOmnimodelPath(selectedNode, availableOmnimodelKeys)
    : null;
  const contextMenuOmnimodel1 = contextMenu ? extractModelFlag(contextMenu.node, 'omnimodel1', LEGACY_MODEL_FIELD_1, LEGACY_MODEL_KEY_1) : false;
  const contextMenuOmnimodel2 = contextMenu ? extractModelFlag(contextMenu.node, 'omnimodel2', LEGACY_MODEL_FIELD_2, LEGACY_MODEL_KEY_2) : false;
  const contextMenuOmnimodelHref = contextMenu && availableOmnimodelKeys
    ? getOmnimodelPath(contextMenu.node, availableOmnimodelKeys)
    : null;

  const nodeByUniqueId = useMemo(() => {
    const map = new Map<string, Node>();
    if (!networkData) return map;
    networkData.nodes.forEach((node) => {
      const uniqueId = getNodeUniqueId(node);
      if (uniqueId) {
        map.set(uniqueId, node);
      }
    });
    return map;
  }, [networkData]);

  const findAncestorNode = useCallback((node: Node | null | undefined): Node | null => {
    if (!node) return null;

    const seen = new Set<string>();
    const uniqueEntries: Array<{ id: string; fallback?: Record<string, unknown> }> = [];

    const pushEntry = (uniqueId?: string, fallback?: Record<string, unknown>) => {
      if (!uniqueId) return;
      if (seen.has(uniqueId)) return;
      seen.add(uniqueId);
      uniqueEntries.push({ id: uniqueId, fallback });
    };

    pushEntry(getNodeUniqueId(node), node.original_model);

    const identicals = node.identical_models;
    if (Array.isArray(identicals)) {
      identicals.forEach((ident) => {
        const uniqueId = (ident as any)?.unique_modelDB_mod_id as string | undefined;
        pushEntry(uniqueId, ident as Record<string, unknown>);
      });
    }

    const candidates: Array<{ node: Node; fallback?: Record<string, unknown> }> = [];
    uniqueEntries.forEach(({ id, fallback }) => {
      const candidate = nodeByUniqueId.get(id);
      if (candidate) {
        candidates.push({ node: candidate, fallback });
      }
    });

    if (candidates.length === 0) {
      return null;
    }

    candidates.sort((a, b) => {
      const icgA = ((a.node.original_model?.ICG ?? (a.fallback as any)?.ICG) ? 0 : 1);
      const icgB = ((b.node.original_model?.ICG ?? (b.fallback as any)?.ICG) ? 0 : 1);
      if (icgA !== icgB) return icgA - icgB;

      const yearA = getCandidateYear(a.node, a.fallback);
      const yearB = getCandidateYear(b.node, b.fallback);
      if (yearA !== yearB) return yearA - yearB;

      const dirA = getCandidateModelDbDir(a.node, a.fallback);
      const dirB = getCandidateModelDbDir(b.node, b.fallback);
      if (dirA !== dirB) return dirA - dirB;

      const uniqueA = getNodeUniqueId(a.node) ?? '';
      const uniqueB = getNodeUniqueId(b.node) ?? '';
      return uniqueA.localeCompare(uniqueB);
    });

    for (const entry of candidates) {
      if (entry.node.id !== node.id) {
        return entry.node;
      }
    }

    return null;
  }, [nodeByUniqueId]);

  const toggleInfoFreeze = useCallback(() => {
    setIsInfoFrozen((prev) => {
      const next = !prev;
      infoFrozenRef.current = next;
      return next;
    });
  }, []);

  useEffect(() => {
    infoFrozenRef.current = isInfoFrozen;
  }, [isInfoFrozen]);

  const clearSelectionsAndDiffs = useCallback(() => {
    commitSelectionState(() => ({ selectedIds: [], sourceIds: [], anchorId: null }));
    setSelectedNode(null);
    setDiffCombinations([]);
    setShowDiffView(false);
    setContextMenu(null);
    setFailedNodes([]);
    setShowFailureBanner(false);
    setIsFailureExpanded(false);
    setIsGeneratingDiffs(false);
    setIsInfoFrozen(false);
    infoFrozenRef.current = false;
    inFlightDiffKeysRef.current.clear();
  }, [commitSelectionState]);

  // Clear only node selections (keep filters and other UI state)
  const clearNodeSelections = useCallback(() => {
    commitSelectionState(() => ({ selectedIds: [], sourceIds: [], anchorId: null }));
    setSelectedNode(null);
  }, [commitSelectionState]);

  // Dark mode removed; keep light theme only

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

  // Close context menu when clicking elsewhere
  useEffect(() => {
    const handleClickOutside = () => setContextMenu(null);
    
    if (contextMenu) {
      document.addEventListener('click', handleClickOutside);
      return () => document.removeEventListener('click', handleClickOutside);
    }
  }, [contextMenu]);

  // Auto-hide banner if no failures remain
  useEffect(() => {
    if (failedNodes.length === 0) {
      setShowFailureBanner(false);
    }
  }, [failedNodes.length]);

  // Retry all failed source fetches
  const retryFailedFetches = async () => {
    if (!networkData) return;
    if (failedNodes.length === 0) return;
    try {
      setRetryingFailed(true);
      const ids = failedNodes.map(f => f.nodeId);
      const { files, failures } = await fetchSourceCode(ids, networkData, fetchedFilesRef.current);
      storeFetchedFiles(files);
      setFailedNodes(failures);
      if (failures.length === 0) setShowFailureBanner(false);
    } finally {
      setRetryingFailed(false);
    }
  };

  

  // Awesome diff generation function
  const generateAwesomeDiffs = useCallback(async (selectedNodes: Node[], useDesignatedSourceTarget: boolean = false) => {
    if (!networkData) {
      return;
    }
    
    setIsGeneratingDiffs(true);
    
    let sourceNodes: Node[] = [];
    let targetNodes: Node[] = [];
    
    if (useDesignatedSourceTarget) {
      sourceNodes = selectionBuckets.sources;
      targetNodes = selectionBuckets.selected;
      if (sourceNodes.length === 0 || targetNodes.length === 0) {
        // Fallback to generic multi-select if designated buckets are empty
        sourceNodes = selectedNodes;
        targetNodes = selectedNodes;
        useDesignatedSourceTarget = false;
      }
    } else {
      sourceNodes = selectedNodes;
      targetNodes = selectedNodes;
    }

    if (sourceNodes.length === 0 && targetNodes.length === 0) {
      alert('Select at least one node to compare.');
      setIsGeneratingDiffs(false);
      return;
    }

    // Enforce selection limit (by unique nodes participating)
    const uniqueNodeIds = new Set<string>([
      ...sourceNodes.map(n => String(n.id)),
      ...targetNodes.map(n => String(n.id)),
    ]);
    if (uniqueNodeIds.size > MAX_SELECTED_NODES) {
      alert('Sorry, please select equal or fewer than 50 nodes.');
      setIsGeneratingDiffs(false);
      return;
    }

    // Quick upper-bound on comparisons to avoid heavy workloads
    const sourcesSet = new Set<string>(sourceNodes.map(n => String(n.id)));
    const targetsSet = new Set<string>(targetNodes.map(n => String(n.id)));
    let intersection = 0;
    for (const id of sourcesSet) if (targetsSet.has(id)) intersection++;
    const estimatedComparisons = (sourcesSet.size * targetsSet.size) - intersection; // skip self-pairs
    if (estimatedComparisons > MAX_COMBINATIONS) {
      alert(`This selection would create ${estimatedComparisons.toLocaleString()} comparisons. Limit each diff run to ${MAX_COMBINATIONS}.`);
      setIsGeneratingDiffs(false);
      return;
    }
    
    const combinationKeys = new Set<string>();
    const combinations: { source: Node; target: Node }[] = [];

    const addCombination = (source: Node, target: Node) => {
      const key = `${String(source.id)}->${String(target.id)}`;
      if (source.id === target.id) return;
      if (combinationKeys.has(key)) return;
      combinationKeys.add(key);
      combinations.push({ source, target });
    };

    generateAllCombinations(sourceNodes, targetNodes).forEach(({ source, target }) => addCombination(source, target));

    const addAncestorPair = (node: Node) => {
      const ancestor = findAncestorNode(node);
      if (ancestor && ancestor.id !== node.id) {
        addCombination(ancestor, node);
      }
    };

    if (useDesignatedSourceTarget) {
      targetNodes.forEach(addAncestorPair);
      sourceNodes.forEach(addAncestorPair);
    } else if (selectedNodes.length === 1) {
      addAncestorPair(selectedNodes[0]);
    }

    if (combinations.length > MAX_COMBINATIONS) {
      alert(`This selection would create ${combinations.length.toLocaleString()} comparisons. Limit each diff run to ${MAX_COMBINATIONS}.`);
      setIsGeneratingDiffs(false);
      return;
    }

    if (combinations.length === 0) {
      alert('Please select at least 2 different nodes, or select a node to trace its ancestors to generate diffs');
      setIsGeneratingDiffs(false);
      return;
    }
    
    // Fetch source code for all involved nodes
    const allNodeIds = Array.from(new Set(
      combinations.flatMap(({ source, target }) => [String(source.id), String(target.id)])
    ));
    const { files: updatedFiles, failures } = await fetchSourceCode(allNodeIds, networkData, fetchedFilesRef.current);
    storeFetchedFiles(updatedFiles);
    setFailedNodes(failures);
    if (failures.length > 0) setShowFailureBanner(true);

    // Seed overlay early so users see progress
    const seeded: { source: Node, target: Node, html?: string }[] = combinations.map(c => ({ ...c }));
    inFlightDiffKeysRef.current.clear();
    setDiffCombinations(seeded);
    setCurrentCombinationIndex(0);
    setShowDiffView(true);
    setContextMenu(null);

    setIsGeneratingDiffs(false);
  }, [networkData, selectionBuckets, storeFetchedFiles, findAncestorNode]);

  // Generate only the comparison being viewed. This keeps large matrices responsive and bounded in memory.
  useEffect(() => {
    if (!showDiffView) return;
    const combo = diffCombinations[currentCombinationIndex];
    if (!combo || combo.html) return;

    const key = `${String(combo.source.id)}->${String(combo.target.id)}`;
    if (inFlightDiffKeysRef.current.has(key)) return;
    inFlightDiffKeysRef.current.add(key);
    setIsGeneratingDiffs(true);

    const sourceCode = fetchedFilesRef.current[String(combo.source.id)]?.source_code;
    const targetCode = fetchedFilesRef.current[String(combo.target.id)]?.source_code;
    const createHtml = async () => {
      if (sourceCode && targetCode) {
        return generateDiffHtml(sourceCode, targetCode, resolvedTheme);
      }
      const missing: string[] = [];
      if (!sourceCode) missing.push(`${combo.source.name || combo.source.id}`);
      if (!targetCode) missing.push(`${combo.target.name || combo.target.id}`);
      return diffPlaceholderDocument(`Missing source code for: ${missing.join(' and ')}`, resolvedTheme);
    };

    void createHtml()
      .then((html) => {
        setDiffCombinations((previous) => {
          const current = previous[currentCombinationIndex];
          if (!current || `${String(current.source.id)}->${String(current.target.id)}` !== key || current.html) {
            return previous;
          }
          const next = [...previous];
          next[currentCombinationIndex] = { ...current, html };
          return next;
        });
      })
      .finally(() => {
        inFlightDiffKeysRef.current.delete(key);
        setIsGeneratingDiffs(false);
      });
  }, [currentCombinationIndex, diffCombinations, resolvedTheme, showDiffView]);

  // Keyboard event handler for shortcuts
  useEffect(() => {
    const handleKeyPress = async (event: KeyboardEvent) => {
      const activeElement = document.activeElement as HTMLElement | null;
      const isEditable = !!activeElement && (
        activeElement.tagName === 'INPUT'
        || activeElement.tagName === 'TEXTAREA'
        || activeElement.tagName === 'SELECT'
        || activeElement.isContentEditable
      );

      if (isEditable && event.code !== 'Escape') return;
      if (sourceView && event.code !== 'Escape') return;
      if (showDiffView && !['ArrowLeft', 'ArrowRight', 'Escape'].includes(event.code)) return;

      if (event.code === 'KeyD') {
        const { sources, selected, all } = selectionBuckets;
        if (sources.length > 0 && selected.length > 0) {
          await generateAwesomeDiffs([...sources, ...selected], true);
        } else if (all.length >= 2) {
          await generateAwesomeDiffs(all, false);
        } else {
          alert('Please select at least 2 nodes to generate diffs');
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
      } else if (event.code === 'KeyA') {
        // Select ancestor subgraph for the focused node
        if (subgraphSelectFn && selectedNode) {
          try { await subgraphSelectFn(selectedNode); } catch {}
        } else {
          alert('Select a node first, then press A to select its ancestor subgraph.');
        }
      } else if (event.code === 'KeyS') {
        if (selectedNode) {
          await openSourceForNode(selectedNode);
        } else {
          alert('Select a node first, then press S to view its source code.');
        }
      } else if (event.code === 'KeyL') {
        toggleInfoFreeze();
      } else if (event.code === 'Escape') {
        // Priority 1: close any modal view
        if (sourceView) {
          closeSourceView();
          return;
        }

        if (showDiffView) {
          setShowDiffView(false);
          setDiffCombinations([]);
          setCurrentCombinationIndex(0);
          setContextMenu(null);
          inFlightDiffKeysRef.current.clear();
          return;
        }

        // Priority 2: cancel any active brush rectangle
        let cancelled = false;
        try {
          if (clearBrushRef.current) {
            clearBrushRef.current();
            cancelled = true;
          }
        } catch {}
        setContextMenu(null);
        if (cancelled) return;

        // Priority 3: if any nodes are selected, clear the selection
        if (selectionBuckets.all.length > 0) {
          clearNodeSelections();
          return;
        }
      } else if (event.code === 'Backspace') {
        event.preventDefault();
        clearSelectionsAndDiffs();
      }
    };

    window.addEventListener('keydown', handleKeyPress);
    return () => window.removeEventListener('keydown', handleKeyPress);
  }, [showDiffView, diffCombinations, fitToViewFunction, subgraphSelectFn, selectedNode, sourceView, openSourceForNode, closeSourceView, toggleInfoFreeze, selectionBuckets, generateAwesomeDiffs, clearNodeSelections, clearSelectionsAndDiffs]);

  useEffect(() => {
    if (availableOmnimodelKeys || (!selectedNode && !contextMenu)) return;
    let cancelled = false;
    async function loadManifest() {
      try {
        const response = await fetch('/api/omnimodels/manifest');
        if (!response.ok) {
          return;
        }
        const data = await response.json();
        const nextKeys = new Set<string>();
        if (data?.channels) {
          for (const channel of data.channels as Array<{ channel: string; models: Array<{ slug: string }> }>) {
            if (!channel?.channel || !Array.isArray(channel.models)) continue;
            for (const model of channel.models) {
              if (!model?.slug) continue;
              nextKeys.add(makeOmnimodelKey({ channel: channel.channel, slug: model.slug }));
            }
          }
        }
        if (!cancelled) {
          setAvailableOmnimodelKeys(nextKeys);
        }
      } catch (error) {
        console.warn('Failed to load omnimodel manifest:', error);
      }
    }
    loadManifest();
    return () => {
      cancelled = true;
    };
  }, [availableOmnimodelKeys, contextMenu, selectedNode]);

  // Load data on mount
  useEffect(() => {
    const controller = new AbortController();

    async function loadData() {
      try {
        const networkResponse = await fetch(INITIAL_NETWORK_DATA_URL, { signal: controller.signal });
        if (!networkResponse.ok) {
          throw new Error(`Network dataset request failed (${networkResponse.status})`);
        }

        const data: unknown = await networkResponse.json();
        setNetworkData(processNetworkDataset(data));
        setDataLoadError(null);
      } catch (error) {
        if (controller.signal.aborted) return;
        console.error('Error loading data:', error);
        setNetworkData(null);
        setDataLoadError('The scientific network dataset could not be loaded. No substitute data is being shown.');
      } finally {
        if (!controller.signal.aborted) {
          setIsLoading(false);
        }
      }
    }
    loadData();
    return () => controller.abort();
  }, []);

  // The initial asset contains every node visible under the default
  // minimum-copies setting. Fetch the complete graph only when the user asks
  // to include one-copy models.
  useEffect(() => {
    if (copiesNumber >= INITIAL_MINIMUM_COPIES || hasFullNetworkData || !networkData) {
      setIsFullNetworkLoading(false);
      return;
    }

    const controller = new AbortController();
    setIsFullNetworkLoading(true);
    setFullNetworkError(null);

    async function loadFullNetwork() {
      try {
        const response = await fetch(FULL_NETWORK_DATA_URL, { signal: controller.signal });
        if (!response.ok) throw new Error(`Network dataset request failed (${response.status})`);
        const data: unknown = await response.json();
        if (controller.signal.aborted) return;
        setNetworkData(processNetworkDataset(data));
        setHasFullNetworkData(true);
      } catch (error) {
        if (controller.signal.aborted) return;
        console.error('Error loading the complete network dataset:', error);
        setCopiesNumber(INITIAL_MINIMUM_COPIES);
        setFullNetworkError('One-copy models could not be loaded. The complete default dataset remains available.');
      } finally {
        if (!controller.signal.aborted) setIsFullNetworkLoading(false);
      }
    }

    loadFullNetwork();
    return () => controller.abort();
  }, [copiesNumber, hasFullNetworkData, networkData]);

  // D3 Force Simulation
  useEffect(() => {
    if (!networkData || !svgRef.current) return;
    if (groupedLayoutCount > 1 && !fixedLocationCircles[groupedLayoutCount]) return;

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

    // Filter data based on current filters (matching original logic)
    const filteredNodes = filterNetworkNodes(
      networkData.nodes,
      selectedIonClasses,
      showICG,
      deferredCopiesNumber,
    );

    const nextVisibleIds = new Set(filteredNodes.map((node) => String(node.id)));
    setVisibleNodeIds((prev) => (areSetsEqual(prev, nextVisibleIds) ? prev : nextVisibleIds));

    commitSelectionState((prev) => {
      const filtered = Array.from(nextVisibleIds);
      const filteredSet = new Set(filtered);
      const selectedIds = prev.selectedIds.filter((id) => filteredSet.has(id));
      const sourceIds = prev.sourceIds.filter((id) => filteredSet.has(id));
      const anchorId = prev.anchorId && filteredSet.has(prev.anchorId) ? prev.anchorId : null;
      return { selectedIds, sourceIds, anchorId };
    });

    const nodeIds = new Set(filteredNodes.map(n => n.id));
    const filteredLinks = networkData.links.filter(link => {
      const sourceId = typeof link.source === 'object' ? link.source.id : link.source;
      const targetId = typeof link.target === 'object' ? link.target.id : link.target;
      const linkWeight = link.weight || link.value || 0;
      // Use > not >= to match original implementation
      return nodeIds.has(sourceId) && nodeIds.has(targetId) && linkWeight > deferredSimilarityScore;
    });

    // Group nodes based on filters
    const uniqueFilterCombinations = getUniqueFilterCombinations(selectedIonClasses, showICG, omnimodel1, omnimodel2);
    const nodesGroupedByFilter = assignNodesToGroups(filteredNodes, uniqueFilterCombinations, selectedIonClasses, showICG, omnimodel1, omnimodel2);
    const nodeIdToLocation = sortAndPositionGroups(nodesGroupedByFilter, fixedLocationCircles);
    
    // Check if we should split groups (more than one group)
    const splitVar = Object.keys(nodesGroupedByFilter).length > 1;
    
    // Create group summaries when groups are split
    if (splitVar) {
      const newGroupSummaries: {key: string, nodeCount: number}[] = [];
      Object.entries(nodesGroupedByFilter).forEach(([groupKey, groupNodes]) => {
        if (groupNodes.length > 0) {
          // Format the group key
          let formattedKey = groupKey;
          
          // Only format if it has multiple parts
          if (groupKey.includes(',')) {
            formattedKey = groupKey
              .replace(',', ': ')
              .replace(/,/g, ' • ')
              .replace(`${LEGACY_MODEL_KEY_1} • ${LEGACY_MODEL_KEY_2}`, 'Omnimodel 1 & 2')
              .replace('Omnimodel 1 • Omnimodel 2', 'Omnimodel 1 & 2');
          }
          
          newGroupSummaries.push({
            key: formattedKey,
            nodeCount: groupNodes.length
          });
        }
      });
      setGroupSummaries(newGroupSummaries);
    } else {
      // Clear summaries when not splitting
      setGroupSummaries([]);
    }

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
    
    // Setup forces exactly as in original simulation.js
    const forceX = splitVar ? 
      d3.forceX((d: any) => nodeIdToLocation[d.id] ? nodeIdToLocation[d.id][0] * width : width / 2).strength(FORCE_STRENGTH) :
      d3.forceX(width / 2);
    
    const forceY = splitVar ?
      d3.forceY((d: any) => nodeIdToLocation[d.id] ? nodeIdToLocation[d.id][1] * height : height / 2).strength(FORCE_STRENGTH) :
      d3.forceY(height / 2);
    
    // Build group label nodes (fixed anchors with collision to keep text readable)
    const labelNodes: any[] = [];
    Object.entries(nodesGroupedByFilter).forEach(([groupKey, groupNodes]) => {
      if (!groupNodes || groupNodes.length === 0) return;
      // Format key similarly to summary
      let formattedKey = groupKey;
      if (groupKey.includes(',')) {
        formattedKey = groupKey
          .replace(',', ': ')
          .replace(/,/g, ' • ')
          .replace(`${LEGACY_MODEL_KEY_1} • ${LEGACY_MODEL_KEY_2}`, 'Omnimodel 1 & 2')
          .replace('Omnimodel 1 • Omnimodel 2', 'Omnimodel 1 & 2');
      }
      // Compute center from positioned node locations (normalized 0..1 from sortAndPositionGroups)
      let sx = 0, sy = 0, c = 0;
      for (const n of groupNodes) {
        const loc = (nodeIdToLocation as any)[n.id];
        if (loc) { sx += loc[0]; sy += loc[1]; c++; }
      }
      // If we don't have explicit group locations (single group case), center label
      const cx = c === 0 ? width / 2 : (sx / c) * width;
      const cy = c === 0 ? height / 2 : (sy / c) * height;
      // Approximate radius from label length to serve as collision obstacle
      const approx = Math.max(26, Math.min(120, (formattedKey.length * 7) / 2 + 18));
      labelNodes.push({ id: `__label__:${groupKey}`, type: 'label', label: formattedKey, x: cx, y: cy, fx: cx, fy: cy, r: approx });
    });

    // Create simulation exactly as in original - note center force is not used in simulation.js
    // but keeping X and Y forces as they handle centering. Include label nodes for collision.
    const simulation = d3.forceSimulation()
      .nodes([...(filteredNodes as any), ...labelNodes])
      .force('charge', chargeForce)
      .force('links', linkForce)
      .force('x', forceX)
      .force('y', forceY)
      .force('collide', d3.forceCollide((d: any) => d.type === 'label' ? d.r : Math.max(2, calculateNodeRadius(d)) + 1).iterations(2));

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
      .attr('r', function (d: Node) {
        const radius = calculateNodeRadius(d);
        d3.select(this).attr('data-base-radius', radius);
        return radius;
      })
      .attr('class', 'graph-node')
      .attr('fill', '#00BFFF')  // DeepSkyBlue - matching original
      .attr('stroke', '#aaa')
      .on('mouseover', function(event: any, d: any) {
        // Update selected node on hover for immediate feedback unless frozen
        if (!infoFrozenRef.current) {
          setSelectedNode(d);
        }
      })
      .on('click', function(event: MouseEvent, d: any) {
        setSelectedNode(d);
        const nodeId = String(d.id);

        commitSelectionState((prev) => {
          const selectedSet = new Set(prev.selectedIds);
          const sourceSet = new Set(prev.sourceIds);
          let anchorId = prev.anchorId;

        if (event.altKey) {
          if (sourceSet.has(nodeId)) {
            sourceSet.delete(nodeId);
          }
          if (selectedSet.has(nodeId)) {
            selectedSet.delete(nodeId);
            if (anchorId === nodeId) {
              const remaining = Array.from(selectedSet);
              const nextAnchor = remaining.length > 0
                ? remaining[remaining.length - 1]
                : (sourceSet.size > 0 ? Array.from(sourceSet)[sourceSet.size - 1] : null);
              anchorId = nextAnchor ?? null;
            }
          }
          return { selectedIds: Array.from(selectedSet), sourceIds: Array.from(sourceSet), anchorId };
        }

        if (event.metaKey || event.ctrlKey) {
          if (selectedSet.has(nodeId)) {
            selectedSet.delete(nodeId);
            if (anchorId === nodeId) {
              const nextAnchor = selectedSet.size > 0 ? Array.from(selectedSet)[selectedSet.size - 1] : (sourceSet.size > 0 ? Array.from(sourceSet)[sourceSet.size - 1] : null);
              anchorId = nextAnchor ?? null;
              }
            } else if (sourceSet.has(nodeId)) {
              sourceSet.delete(nodeId);
              selectedSet.add(nodeId);
              anchorId = nodeId;
            } else {
              selectedSet.add(nodeId);
              anchorId = nodeId;
            }
            return { selectedIds: Array.from(selectedSet), sourceIds: Array.from(sourceSet), anchorId };
          }

        if (event.shiftKey) {
          if (sourceSet.has(nodeId)) {
            sourceSet.delete(nodeId);
          }
          if (selectedSet.has(nodeId)) {
            selectedSet.delete(nodeId);
            if (anchorId === nodeId) {
              const remaining = Array.from(selectedSet);
              const nextAnchor = remaining.length > 0
                ? remaining[remaining.length - 1]
                : (sourceSet.size > 0 ? Array.from(sourceSet)[sourceSet.size - 1] : null);
              anchorId = nextAnchor ?? null;
            }
          } else {
            selectedSet.add(nodeId);
            if (!anchorId) {
              anchorId = nodeId;
            }
          }
          return { selectedIds: Array.from(selectedSet), sourceIds: Array.from(sourceSet), anchorId };
        }

          return { selectedIds: [nodeId], sourceIds: [], anchorId: nodeId };
        });
      })
      .on('contextmenu', function(event: any, d: any) {
        event.preventDefault();
        
        const x = event.clientX;
        const y = event.clientY;
        
        setContextMenu({ x, y, node: d });
      })
      .call(d3.drag<any, any>()
        .on('start', dragstarted)
        .on('drag', dragged)
        .on('end', dragended) as any);

    const collectIdsInRect = (x0: number, y0: number, x1: number, y1: number): string[] => {
      const ids: string[] = [];
      node.each(function (d: any) {
        const nx = typeof d.x === 'number' ? d.x : 0;
        const ny = typeof d.y === 'number' ? d.y : 0;
        if (x0 <= nx && nx <= x1 && y0 <= ny && ny <= y1) {
          ids.push(String(d.id));
        }
      });
      return uniqueIdList(ids);
    };

    const applyBrushSelection = (mode: 'add' | 'subtract' | 'replace', ids: string[]) => {
      if (ids.length === 0 && mode === 'replace') {
        commitSelectionState(() => ({ selectedIds: [], sourceIds: [], anchorId: null }));
        return;
      }

      commitSelectionState((prev) => {
        const selectedSet = new Set(prev.selectedIds);
        const sourceSet = new Set(prev.sourceIds);

        if (mode === 'subtract') {
          ids.forEach((id) => {
            if (!sourceSet.has(id)) {
              selectedSet.delete(id);
            }
          });
          return { selectedIds: Array.from(selectedSet), sourceIds: Array.from(sourceSet), anchorId: prev.anchorId };
        }

        if (mode === 'add') {
          ids.forEach((id) => {
            if (sourceSet.has(id)) {
              sourceSet.delete(id);
            }
            selectedSet.add(id);
          });
          const anchorId = prev.anchorId ?? (ids[ids.length - 1] ?? null);
          return { selectedIds: Array.from(selectedSet), sourceIds: Array.from(sourceSet), anchorId };
        }

        const anchorId = ids[ids.length - 1] ?? null;
        return { selectedIds: ids, sourceIds: [], anchorId };
      });
    };

    const marqueeLayer = g.append('g').attr('class', 'selection-marquee-layer').style('pointer-events', 'none');
    const marqueeRect = marqueeLayer
      .append('rect')
      .attr('class', 'selection-marquee')
      .attr('fill', 'rgba(59,130,246,0.18)')
      .attr('stroke', '#2563eb')
      .attr('stroke-width', 1)
      .attr('rx', 6)
      .attr('ry', 6)
      .attr('visibility', 'hidden');
    marqueeRectRef.current = marqueeRect;

    const computeMode = (): 'replace' | 'add' | 'subtract' => {
      const { alt, meta, ctrl, shift } = modifierStateRef.current;
      if (alt) return 'subtract';
      if (meta || ctrl || shift) return 'add';
      return 'replace';
    };

    const svgElement = svg.node();
    const gElement = g.node();

    const handlePointerDown = (event: PointerEvent) => {
      if (event.button !== 0 || !gElement) return;
      const target = event.target as Element | null;
      const mode = computeMode();
      const onNode = !!target?.closest('circle.graph-node');
      if (onNode && mode === 'replace') {
        return;
      }

      event.preventDefault();
      const [startX, startY] = d3.pointer(event, gElement);
      marqueeStateRef.current = { active: true, startX, startY, mode };
      isBrushingRef.current = true;
      marqueeRect
        .attr('visibility', 'visible')
        .attr('x', startX)
        .attr('y', startY)
        .attr('width', 0)
        .attr('height', 0);

      const handlePointerMove = (moveEvent: PointerEvent) => {
        if (!marqueeStateRef.current.active) return;
        marqueeStateRef.current.mode = computeMode();
        const [mx, my] = d3.pointer(moveEvent, gElement);
        const { startX: sx, startY: sy } = marqueeStateRef.current;
        marqueeRect
          .attr('x', Math.min(sx, mx))
          .attr('y', Math.min(sy, my))
          .attr('width', Math.abs(mx - sx))
          .attr('height', Math.abs(my - sy));
      };

      const finalizeSelection = (upEvent: PointerEvent) => {
        window.removeEventListener('pointermove', handlePointerMove);
        window.removeEventListener('pointerup', finalizeSelection);
        if (!marqueeStateRef.current.active) return;

        const { startX: sx, startY: sy, mode: finalMode } = marqueeStateRef.current;
        marqueeStateRef.current.active = false;
        marqueeRect.attr('visibility', 'hidden');
        isBrushingRef.current = false;

        const [ex, ey] = d3.pointer(upEvent, gElement);
        const x0 = Math.min(sx, ex);
        const y0 = Math.min(sy, ey);
        const x1 = Math.max(sx, ex);
        const y1 = Math.max(sy, ey);
        const width = Math.abs(ex - sx);
        const height = Math.abs(ey - sy);

        if (width < 2 && height < 2) {
          clearBrushRef.current = null;
          return;
        }

        const ids = collectIdsInRect(x0, y0, x1, y1);
        applyBrushSelection(finalMode, ids);
        setSelectedNode(null);
        clearBrushRef.current = null;
      };

      window.addEventListener('pointermove', handlePointerMove);
      window.addEventListener('pointerup', finalizeSelection);

      clearBrushRef.current = () => {
        marqueeStateRef.current.active = false;
        marqueeRect.attr('visibility', 'hidden');
        isBrushingRef.current = false;
        window.removeEventListener('pointermove', handlePointerMove);
        window.removeEventListener('pointerup', finalizeSelection);
        clearBrushRef.current = null;
      };
    };

    if (svgElement) {
      svgElement.addEventListener('pointerdown', handlePointerDown, { capture: true });
    }

    // Group label visuals (text + invisible circle for reference)
    const labelsG = g.append('g').attr('class', 'group-labels');
    const labelSel = labelsG.selectAll('g')
      .data(labelNodes)
      .enter()
      .append('g')
      .attr('class', 'group-label');
    // Optional invisible circle for debugging/collision reference
    labelSel.append('circle')
      .attr('r', (d: any) => d.r)
      .attr('fill', 'white')
      .attr('fill-opacity', 0)
      .attr('stroke', 'none')
      .attr('pointer-events', 'none');
    // Text label with halo for readability
    labelSel.append('text')
      .attr('text-anchor', 'middle')
      .attr('dominant-baseline', 'middle')
      .attr('font-size', 12)
      .attr('font-weight', 600)
      .attr('stroke', '#ffffff')
      .attr('stroke-width', 3)
      .attr('paint-order', 'stroke')
      .attr('fill', '#0f172a')
      .text((d: any) => d.label)
      .attr('pointer-events', 'none');

    // Add tooltips
    node.append('title')
      .text((d: Node) => `${getNodeDisplayName(d)}\nClass: ${formatIonClassForDisplay(d.original_model?.ion_class || d.ion_class)}\nICG: ${d.icg ? 'Yes' : 'No'}`);

    syncSelectionToDom(selectionStateRef.current);

    // Calculate bounding box of all nodes for proper zoom setup
    let cachedBounds: { minX: number; maxX: number; minY: number; maxY: number } | null = null;
    const calculateBounds = () => {
      if (cachedBounds) return cachedBounds;
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
      
      cachedBounds = { minX, maxX, minY, maxY };
      return cachedBounds;
    };

    // Calculate dynamic zoom constraints
    const calculateMinScale = () => {
      const bounds = calculateBounds();
      const boundsWidth = bounds.maxX - bounds.minX;
      const boundsHeight = bounds.maxY - bounds.minY;
      
      // Calculate scale that fits all nodes with margin
      return Math.min(
        width / boundsWidth * 0.85,  // 85% to ensure all nodes visible
        height / boundsHeight * 0.85
      );
    };
    
    // Enhanced zoom behavior with progressive centering
    const zoom = d3.zoom()
      .filter((event: any) => {
        const e: any = event;
        const t = e?.type as string | undefined;
        // Always allow wheel zoom
        if (t && (t === 'wheel' || t === 'mousewheel')) return true;
        // If interacting with the brush overlay/selection, block zoom
        try {
          const target = e?.target as Element | null | undefined;
          if (target && target.closest('g.brush')) return false;
        } catch {}
        // If a brush is (about to be) used via modifiers, suppress pan-zoom
        if (e && (e.shiftKey || e.ctrlKey || e.metaKey || e.altKey)) return false;
        // Block pan/zoom while actively brushing
        if (isBrushingRef.current) return false;
        return true;
      })
      .scaleExtent([0.01, 10])  // Start with very permissive range
      .on('zoom', (event) => {
        // If actively brushing, ignore any pan/zoom updates
        if (isBrushingRef.current) return;
        let { x, y, k } = event.transform;
        const bounds = calculateBounds();
        
        // Calculate center of bounds (center of gravity)
        const centerX = (bounds.minX + bounds.maxX) / 2;
        const centerY = (bounds.minY + bounds.maxY) / 2;
        
        // Calculate minimum scale dynamically
        const minScale = calculateMinScale();
        
        // Progressive centering as zoom approaches minimum
        if (k <= minScale * 1.5) {
          // Calculate ideal centered position
          const idealX = width / 2 - centerX * k;
          const idealY = height / 2 - centerY * k;
          
          if (k <= minScale) {
            // Hard constraint at minimum zoom
            k = minScale;
            x = idealX;
            y = idealY;
          } else {
            // Progressive centering between minScale and minScale * 1.5
            // The closer to minScale, the stronger the centering force
            const centeringFactor = 1 - ((k - minScale) / (minScale * 0.5));
            const smoothFactor = Math.pow(centeringFactor, 2); // Smooth curve
            
            // Blend current position with ideal centered position
            x = x * (1 - smoothFactor) + idealX * smoothFactor;
            y = y * (1 - smoothFactor) + idealY * smoothFactor;
          }
          
          const constrainedTransform = d3.zoomIdentity
            .translate(x, y)
            .scale(k);
          
          g.attr('transform', constrainedTransform.toString());
          svg.property('__zoom', constrainedTransform);
        } else {
          // Normal zoom behavior when well above minimum
          g.attr('transform', event.transform);
        }
      });

    svg.call(zoom as any);
    

    // Function to fit visualization to view
    const fitToView = () => {
      const bounds = calculateBounds();
      const boundsWidth = bounds.maxX - bounds.minX;
      const boundsHeight = bounds.maxY - bounds.minY;
      
      // Use same calculation as minimum scale for consistency
      const scale = Math.min(
        width / boundsWidth * 0.85,
        height / boundsHeight * 0.85
      );
      
      // Calculate center of bounds (center of gravity)
      const centerX = (bounds.minX + bounds.maxX) / 2;
      const centerY = (bounds.minY + bounds.maxY) / 2;
      
      // Calculate translation to center the content
      const translateX = width / 2 - centerX * scale;
      const translateY = height / 2 - centerY * scale;
      
      // Apply the transform
      const transform = d3.zoomIdentity
        .translate(translateX, translateY)
        .scale(scale);
        
      svg.transition()
        .duration(750)
        .call(zoom.transform as any, transform);
    };

    // Store the fit function for keyboard access
    setFitToViewFunction(() => fitToView);

    // Add fit-to-view functionality on double-click
    svg.on('dblclick.zoom', null); // Remove default double-click zoom
    svg.on('dblclick', fitToView);

    // Provide subgraph selection function for keyboard 'A'/'S'
    const selectAncestorSubgraph = async (d: any) => {
      // Build subgraph from provided node
      let subgraphNodeIds = [d.id];
      let previousSize = -1;
      let currentSize = 1;
      while (previousSize !== currentSize) {
        previousSize = currentSize;
        subgraphNodeIds = extendSubgraphByConnectedNodes(filteredLinks, subgraphNodeIds);
        currentSize = subgraphNodeIds.length;
      }

      // Nodes in subgraph
      const nodesInSubgraph = node.filter((nodeData: any) => subgraphNodeIds.includes(nodeData.id));
      const subgraphData = nodesInSubgraph.data();
      const yearsInSubgraph = subgraphData.map((nodeData: any) => nodeData.original_model?.Year || 2023);
      const earliestYear = Math.min(...yearsInSubgraph);
      let originalNodes = subgraphData.filter((nodeData: any) => (nodeData.original_model?.Year || 2023) === earliestYear);
      if (originalNodes.length > 1) {
        const modelDBDirs = originalNodes.map((nodeData: any) => parseInt(nodeData.original_model?.modelDB_dir || '999999'));
        const smallestModelDBDir = Math.min(...modelDBDirs);
        originalNodes = originalNodes.filter((nodeData: any) => parseInt(nodeData.original_model?.modelDB_dir || '999999') === smallestModelDBDir);
      }
      const sourceIds = originalNodes.map((nodeData: any) => nodeData.id);
      const newSourceIds = sourceIds.map((id) => String(id));
      const newTargetIds = subgraphNodeIds
        .filter((id) => !sourceIds.includes(id))
        .map((id) => String(id));

      commitSelectionState(() => ({
        selectedIds: newTargetIds,
        sourceIds: newSourceIds,
        anchorId: newTargetIds[newTargetIds.length - 1] ?? newSourceIds[newSourceIds.length - 1] ?? null,
      }));
      setSelectedNode(null);

      // Fetch source code for all nodes in subgraph
      if (networkData) {
        const allNodeIds = [...newSourceIds, ...newTargetIds];
        const { files: updatedFetchedFiles, failures } = await fetchSourceCode(allNodeIds, networkData, fetchedFilesRef.current);
        storeFetchedFiles(updatedFetchedFiles);
        setFailedNodes(failures);
        if (failures.length > 0) setShowFailureBanner(true);
      }

    };
    setSubgraphSelectFn(() => selectAncestorSubgraph);

    // Simulation tick (matching original - no boundary constraints)
    let tickFrame: number | null = null;
    simulation.on('tick', () => {
      if (tickFrame !== null) return;
      tickFrame = window.requestAnimationFrame(() => {
        tickFrame = null;
        link
          .attr('x1', (d: any) => d.source.x)
          .attr('y1', (d: any) => d.source.y)
          .attr('x2', (d: any) => d.target.x)
          .attr('y2', (d: any) => d.target.y);

        node
          .attr('cx', (d: any) => d.x)
          .attr('cy', (d: any) => d.y);
        g.selectAll<SVGCircleElement, any>('g.glow-layer circle.glow-halo')
          .attr('cx', (d: any) => d.x)
          .attr('cy', (d: any) => d.y);
        labelSel.attr('transform', (d: any) => `translate(${d.x},${d.y})`);
      });
    });
    simulation.on('end.bounds', () => {
      cachedBounds = null;
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
      if (tickFrame !== null) {
        window.cancelAnimationFrame(tickFrame);
      }
      if (clearBrushRef.current) {
        try { clearBrushRef.current(); } catch {}
        clearBrushRef.current = null;
      }
      marqueeRectRef.current = null;
      if (svgElement) {
        svgElement.removeEventListener('pointerdown', handlePointerDown, { capture: true } as any);
      }
    };
  }, [networkData, fixedLocationCircles, groupedLayoutCount, selectedIonClasses, deferredSimilarityScore, showICG, omnimodel1, omnimodel2, deferredCopiesNumber, containerDimensions, resolvedTheme, commitSelectionState, syncSelectionToDom, storeFetchedFiles]);

  useEffect(() => {
    if (!networkData) {
      setSearchMatches((prev) => (prev.size === 0 ? prev : new Set()));
      return;
    }

    const term = normalizeQueryTerm(searchTerm);
    if (term.length < 3) {
      setSearchMatches((prev) => (prev.size === 0 ? prev : new Set()));
      return;
    }

    const visibleIds = visibleNodeIds;
    const restrictToVisible = visibleIds.size > 0;
    const matches = new Set<string>();

    for (const node of networkData.nodes) {
      const nodeId = String(node.id);
      if (restrictToVisible && !visibleIds.has(nodeId)) continue;

      const displayName = getNodeDisplayName(node);
      const uniqueId = getNodeUniqueId(node);
      const label = typeof node.label === 'string' ? node.label : null;
      const nodeName = typeof node.name === 'string' ? node.name : null;
      const modFilename = typeof node.original_model?.mod_filename === 'string' ? node.original_model.mod_filename : null;
      const modFilepath = typeof node.original_model?.mod_filepath === 'string' ? node.original_model.mod_filepath : null;
      const identicalEntries = Array.isArray(node.identical_models) ? node.identical_models : [];

      const haystacks = [
        nodeId,
        label,
        nodeName,
        displayName,
        uniqueId,
        modFilename,
        modFilepath,
        ...identicalEntries.flatMap((entry) => {
          if (!entry || typeof entry !== 'object') return [];
          const castEntry = entry as Record<string, unknown>;
          const parts: string[] = [];
          const id = castEntry.unique_modelDB_mod_id;
          const filename = castEntry.mod_filename;
          const filepath = castEntry.mod_filepath;
          if (typeof id === 'string') parts.push(id);
          if (typeof filename === 'string') parts.push(filename);
          if (typeof filepath === 'string') parts.push(filepath);
          return parts;
        }),
      ].filter((value): value is string => Boolean(value))
       .map((value) => value.toLowerCase());

      if (haystacks.some((field) => field.includes(term))) {
        matches.add(nodeId);
      }
    }

    setSearchMatches((prev) => (areSetsEqual(prev, matches) ? prev : matches));
  }, [searchTerm, networkData, visibleNodeIds]);

  useEffect(() => {
    if (!svgRef.current) return;

    const svg = d3.select(svgRef.current);
    const nodes = svg.selectAll<SVGCircleElement, Node>('.graph-node');
    if (nodes.empty()) return;

    const highlightedIds = searchMatches;
    const baseFill = resolvedTheme === 'dark' ? '#3b82f6' : '#00BFFF';
    const baseStroke = resolvedTheme === 'dark' ? '#475569' : '#aaa';

    // Build map of base radii for nodes
    const radii = new Map<string, number>();
    nodes.each(function (d) {
      const circle = d3.select<SVGCircleElement, Node>(this);
      const baseRadius = Number(circle.attr('data-base-radius')) || Number(circle.attr('r')) || 4;
      radii.set(String(d.id), baseRadius);
    });

    // Ensure a glow layer exists behind nodes (between links and nodes)
    const gRoot = svg.select('g');
    let glowLayer = gRoot.select<SVGGElement>('g.glow-layer');
    if (glowLayer.empty()) {
      glowLayer = gRoot.append('g').attr('class', 'glow-layer');
      // Place it beneath nodes group if present
      const nodesGroup = gRoot.select('g.nodes');
      if (!nodesGroup.empty()) {
        // Move glowLayer just before nodes group
        // D3 doesn't have a direct insert for existing nodes; recreate order by lowering nodes
        glowLayer.lower();
      }
    }

    // Collect matched node data for halos
    const matchedData: Node[] = [];
    nodes.each(function (d) {
      if (highlightedIds.has(String(d.id))) matchedData.push(d);
    });

    const halos = glowLayer.selectAll<SVGCircleElement, Node>('circle.glow-halo').data(matchedData, (d: any) => String(d.id));
    halos.exit().remove();
    const halosEnter = halos.enter().append('circle').attr('class', 'glow-halo search-breathe');
    const halosMerged = halosEnter.merge(halos);
    halosMerged
      .attr('cx', (d: any) => (typeof d.x === 'number' ? d.x : 0))
      .attr('cy', (d: any) => (typeof d.y === 'number' ? d.y : 0))
      .attr('r', (d: any) => {
        const baseRadius = radii.get(String(d.id)) || 4;
        const extra = Math.max(5, Math.min(14, baseRadius * 0.6));
        return baseRadius + extra;
      })
      .attr('fill', 'rgba(34, 197, 94, 0.35)')
      .attr('stroke', 'transparent');

    // Style node circles; apply green inner color for matches (except source nodes)
    nodes.each(function (d) {
      const circle = d3.select<SVGCircleElement, Node>(this);
      const nodeId = String(d.id);
      const isMatch = highlightedIds.has(nodeId);
      const isSource = circle.classed('source-node');
      const isSelected = circle.classed('selected-node') || isSource;
      const baseWidth = isSelected ? 2 : 1;
      const baseRadius = radii.get(nodeId) || Number(circle.attr('r')) || 4;

      if (isMatch) {
        circle
          .classed('search-hit', true)
          .attr('stroke-width', baseWidth)
          .attr('r', baseRadius)
          .style('filter', null);
        // Override inner color to green except for source nodes (keep yellow for sources)
        circle.classed('search-override', !isSource);
        circle.raise();
      } else {
        circle
          .classed('search-hit', false)
          .attr('stroke-width', baseWidth)
          .attr('r', baseRadius)
          .style('filter', null);
        circle.classed('search-override', false);
        if (!isSelected) {
          circle.attr('fill', baseFill).attr('stroke', baseStroke);
        }
      }
    });
  }, [searchMatches, resolvedTheme, visibleNodeIds, containerDimensions]);

  if (isLoading) {
    return (
      <div className="flex items-center justify-center h-screen bg-slate-50 dark:bg-slate-900 text-slate-900 dark:text-slate-100 px-6">
        <div className="w-64 text-center space-y-4">
          <div className="text-lg font-semibold">Loading Visualizer</div>
          <div className="w-full h-1.5 bg-slate-200 dark:bg-slate-700 rounded-full overflow-hidden">
            <div className="h-full w-full bg-slate-500 dark:bg-slate-400 rounded-full origin-left animate-[indeterminate_1.5s_ease-in-out_infinite]" />
          </div>
          <style>{`@keyframes indeterminate { 0% { transform: scaleX(0); transform-origin: left; } 50% { transform: scaleX(1); transform-origin: left; } 50.1% { transform: scaleX(1); transform-origin: right; } 100% { transform: scaleX(0); transform-origin: right; } }`}</style>
        </div>
      </div>
    );
  }

  if (dataLoadError || !networkData) {
    return (
      <div className="flex items-center justify-center h-screen bg-slate-50 dark:bg-slate-900 p-8">
        <div role="alert" className="text-center max-w-lg rounded-xl border border-red-200 dark:border-red-900 bg-white dark:bg-slate-800 p-8 shadow-sm">
          <h1 className="text-2xl font-bold text-slate-900 dark:text-white mb-4">
            Visualizer unavailable
          </h1>
          <p className="text-slate-700 dark:text-slate-300 mb-6">
            {dataLoadError || 'The scientific network dataset is unavailable. No substitute data is being shown.'}
          </p>
          <div className="flex items-center justify-center gap-3">
            <button
              type="button"
              onClick={() => window.location.reload()}
              className="px-5 py-2.5 rounded-lg bg-blue-600 text-white hover:bg-blue-700 transition-colors"
            >
              Retry
            </button>
            <NextLink
              href="/"
              className="px-5 py-2.5 rounded-lg border border-slate-300 dark:border-slate-600 text-slate-800 dark:text-slate-100 hover:bg-slate-100 dark:hover:bg-slate-700 transition-colors"
            >
              Return home
            </NextLink>
          </div>
        </div>
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
            The ICG Explorer requires a larger screen for the best experience.
            Please access this tool on a tablet, laptop, or desktop computer.
          </p>
          <NextLink
            href="/"
            className="inline-block px-6 py-3 rounded-lg bg-blue-600 text-white hover:bg-blue-700 transition-colors"
          >
            Return to Home
          </NextLink>
        </div>
      </div>
    );
  }

  return (
    <div className="h-screen bg-slate-50 dark:bg-slate-900 overflow-hidden scrollbar-none">
      {/* Header */}
      <div className="border-b border-slate-200 dark:border-slate-700 bg-white dark:bg-slate-800">
        <div className="px-4 py-3">
          <div className="flex flex-wrap items-center justify-between gap-3 md:gap-4">
            <h1 className="text-xl font-bold text-slate-900 dark:text-white">
              ICG Explorer
            </h1>
            <form
              onSubmit={(event) => event.preventDefault()}
              className="order-3 w-full md:order-2 md:w-auto md:flex-1 md:max-w-lg"
            >
              <label htmlFor="network-search" className="sr-only">
                Search nodes
              </label>
              <div className="relative">
                <span className="pointer-events-none absolute inset-y-0 left-3 flex items-center text-slate-400 dark:text-slate-500">
                  <svg className="h-4 w-4" viewBox="0 0 24 24" fill="none" stroke="currentColor">
                    <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2} d="M21 21l-4.35-4.35m1.35-4.65a6 6 0 11-12 0 6 6 0 0112 0z" />
                  </svg>
                </span>
                <input
                  id="network-search"
                  type="search"
                  value={searchTerm}
                  onChange={(event) => setSearchTerm(event.target.value)}
                  placeholder="Search nodes, IDs, or files..."
                  className="w-full rounded-md border border-slate-300 bg-white py-2 pl-10 pr-24 text-sm text-slate-900 placeholder-slate-400 shadow-sm focus:border-blue-500 focus:outline-none focus:ring-2 focus:ring-blue-500 dark:border-slate-600 dark:bg-slate-900 dark:text-slate-100 dark:placeholder-slate-400 dark:focus:border-blue-400 dark:focus:ring-blue-400"
                />
                {searchTerm.length > 0 && (
                  <div className="absolute inset-y-0 right-2 flex items-center gap-2">
                    <span className="pointer-events-none select-none rounded-full bg-slate-100 px-2 py-0.5 text-xs font-medium text-slate-600 dark:bg-slate-800 dark:text-slate-300">
                      {searchMatches.size === 0
                        ? 'No hits'
                        : searchMatches.size === 1
                          ? '1 hit'
                          : `${searchMatches.size} hits`}
                    </span>
                    <button
                      type="button"
                      onClick={() => setSearchTerm('')}
                      className="flex h-6 w-6 items-center justify-center rounded-full text-slate-400 transition-colors hover:bg-slate-200 hover:text-slate-600 dark:text-slate-500 dark:hover:bg-slate-700 dark:hover:text-slate-300"
                      aria-label="Clear search"
                    >
                      <svg className="h-3.5 w-3.5" viewBox="0 0 24 24" fill="none" stroke="currentColor">
                        <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2} d="M6 18L18 6M6 6l12 12" />
                      </svg>
                    </button>
                  </div>
                )}
              </div>
            </form>
            <div className="order-2 flex flex-wrap items-center gap-2 md:order-3 md:gap-3">
              <NextLink
                href="/"
                className="px-4 py-2 text-sm rounded-lg border border-slate-300 dark:border-slate-600 hover:bg-slate-100 dark:hover:bg-slate-700 transition-colors cursor-pointer"
              >
                Back to Home
              </NextLink>
              {/* Removed nav actions: Copy copies, Export CSV */}
            </div>
          </div>
        </div>
      </div>

      {/* Failed fetch banner */}
      {showFailureBanner && failedNodes.length > 0 && (
        <div className="fixed top-[57px] left-0 right-0 z-40">
          <div className="mx-4 mt-3 bg-amber-50 dark:bg-amber-900/30 border border-amber-200 dark:border-amber-800 rounded-lg shadow p-3">
            <div className="flex items-center justify-between">
              <div className="flex items-center space-x-2 text-amber-900 dark:text-amber-100 text-sm">
                <svg className="w-4 h-4" viewBox="0 0 24 24" fill="none" stroke="currentColor"><path strokeLinecap="round" strokeLinejoin="round" strokeWidth="2" d="M12 9v4m0 4h.01M10.29 3.86L1.82 18a2 2 0 001.71 3h16.94a2 2 0 001.71-3L13.71 3.86a2 2 0 00-3.42 0z"/></svg>
                <span>Failed to load {failedNodes.length} source {failedNodes.length === 1 ? 'file' : 'files'}.</span>
              </div>
              <div className="flex items-center space-x-2">
                <button
                  onClick={() => setIsFailureExpanded(v => !v)}
                  className="px-2 py-1 text-xs rounded bg-amber-100 dark:bg-amber-800 text-amber-900 dark:text-amber-100 hover:bg-amber-200 dark:hover:bg-amber-700 cursor-pointer"
                >
                  {isFailureExpanded ? 'Hide details' : 'Show details'}
                </button>
                <button
                  onClick={retryFailedFetches}
                  disabled={retryingFailed}
                  className={`px-2 py-1 text-xs rounded cursor-pointer ${retryingFailed ? 'opacity-60 cursor-not-allowed' : ''} bg-blue-600 text-white hover:bg-blue-700`}
                >
                  {retryingFailed ? 'Retrying…' : 'Retry all'}
                </button>
                <button
                  onClick={() => setShowFailureBanner(false)}
                  className="px-2 py-1 text-xs rounded bg-slate-200 dark:bg-slate-700 text-slate-800 dark:text-slate-100 hover:bg-slate-300 dark:hover:bg-slate-600 cursor-pointer"
                >
                  Dismiss
                </button>
              </div>
            </div>
            {isFailureExpanded && (
              <div className="mt-2 max-h-40 overflow-auto scrollbar-none text-sm">
                <ul className="list-disc pl-5 text-amber-900 dark:text-amber-100">
                  {failedNodes.map((f, idx) => (
                    <li key={`${f.nodeId}-${idx}`} className="flex items-center justify-between pr-2">
                      <span>
                        {f.name || f.nodeId}
                        {f.uniqueId ? ` (${f.uniqueId})` : ''}
                        {f.error ? ` – ${f.error}` : ''}
                      </span>
                      <button
                        className="text-xs underline hover:no-underline"
                        onClick={async () => {
                          if (!networkData) return;
                          setRetryingFailed(true);
                          try {
                            const { files, failures } = await fetchSourceCode([f.nodeId], networkData, fetchedFilesRef.current);
                            storeFetchedFiles(files);
                            const stillFailed = failures.find(ff => ff.nodeId === f.nodeId);
                            setFailedNodes(prev => {
                              const others = prev.filter(p => p.nodeId !== f.nodeId);
                              return stillFailed ? [...others, stillFailed] : others;
                            });
                          } finally {
                            setRetryingFailed(false);
                          }
                        }}
                      >Retry</button>
                    </li>
                  ))}
                </ul>
              </div>
            )}
          </div>
        </div>
      )}

      <div className="flex h-[calc(100vh-57px)]">
        {/* Left Control Panel */}
        <div className="w-80 bg-white dark:bg-slate-800 border-r border-slate-200 dark:border-slate-700 overflow-y-auto scrollbar-none flex-shrink-0 viz-left-panel" style={{ ['--viz-scale' as any]: leftFontScale }}>
          <div className="p-4 space-y-6">
            {/* Ion Channel Class Filter */}
            <div>
              <h3 id="ion-class-filter-label" className="text-sm font-semibold text-slate-700 dark:text-slate-300 mb-3">
                Ion Channel Class
                {!selectedIonClasses.has('all') && selectedIonClasses.size > 0 && (
                  <span className="ml-2 text-xs font-normal text-slate-500 dark:text-slate-400">
                    ({selectedIonClasses.size} selected)
                  </span>
                )}
              </h3>
              <div className="grid grid-cols-3 gap-2" role="group" aria-labelledby="ion-class-filter-label">
                {ION_CLASS_FILTER_OPTIONS.map((cls) => (
                  <button
                    type="button"
                    key={cls}
                    aria-pressed={selectedIonClasses.has(cls)}
                    onClick={() => {
                      const newSelection = new Set(selectedIonClasses);
                      if (cls === 'all') {
                        // If 'all' is clicked, clear selection and select only 'all'
                        setSelectedIonClasses(new Set(['all']));
                      } else {
                        // Toggle the specific class
                        if (newSelection.has(cls)) {
                          newSelection.delete(cls);
                          // If nothing is selected, default to 'all'
                          if (newSelection.size === 0 || (newSelection.size === 1 && newSelection.has('all'))) {
                            setSelectedIonClasses(new Set(['all']));
                          }
                        } else {
                          // Remove 'all' when selecting specific classes
                          newSelection.delete('all');
                          newSelection.add(cls);
                        }
                        setSelectedIonClasses(newSelection);
                      }
                    }}
                    className={`px-3 py-1.5 text-sm rounded-md transition-colors cursor-pointer ${
                      selectedIonClasses.has(cls)
                        ? 'bg-blue-600 text-white'
                        : 'bg-slate-100 dark:bg-slate-700 text-slate-700 dark:text-slate-300 hover:bg-slate-200 dark:hover:bg-slate-600'
                    }`}
                  >
                    {cls === 'all' ? 'All' : cls}
                  </button>
                ))}
              </div>
            </div>

            {/* Copies Filter */}
            <div>
              <h3 className="text-sm font-semibold text-slate-700 dark:text-slate-300 mb-3">
                Minimum Copies: {copiesNumber}
              </h3>
              <div className="flex items-center space-x-3">
                <button
                  type="button"
                  onClick={() => setCopiesNumber(Math.max(1, copiesNumber - 1))}
                  aria-label="Decrease minimum copies"
                  className="w-8 h-8 rounded-md bg-slate-100 dark:bg-slate-700 hover:bg-slate-200 dark:hover:bg-slate-600 flex items-center justify-center cursor-pointer"
                >
                  −
                </button>
                <input
                  id="minimum-copies"
                  type="range"
                  min="1"
                  max="10"
                  value={copiesNumber}
                  onChange={(e) => setCopiesNumber(Number(e.target.value))}
                  aria-label="Minimum copies"
                  aria-valuetext={`${copiesNumber} copies`}
                  className="flex-1 cursor-pointer"
                />
                <button
                  type="button"
                  onClick={() => setCopiesNumber(Math.min(10, copiesNumber + 1))}
                  aria-label="Increase minimum copies"
                  className="w-8 h-8 rounded-md bg-slate-100 dark:bg-slate-700 hover:bg-slate-200 dark:hover:bg-slate-600 flex items-center justify-center cursor-pointer"
                >
                  +
                </button>
                <span className="text-sm text-slate-600 dark:text-slate-400 min-w-[20px] text-center">
                  {copiesNumber}
                </span>
              </div>
              {isFullNetworkLoading && (
                <p role="status" className="mt-2 text-xs text-slate-500 dark:text-slate-400">
                  Loading one-copy models…
                </p>
              )}
              {fullNetworkError && (
                <p role="alert" className="mt-2 text-xs text-red-700 dark:text-red-400">
                  {fullNetworkError}
                </p>
              )}
            </div>

            {/* Minimum code similarity (persisted in ?sim=) */}
            <div>
              <h3 className="text-sm font-semibold text-slate-700 dark:text-slate-300 mb-3">
                Minimum code similarity <span className="ml-2 inline-flex items-center rounded-full border border-slate-300 px-2 py-0.5 text-xs text-slate-700 dark:text-slate-300">≥ {similarityScore}%</span>
              </h3>
              <div className="flex items-center space-x-3">
                <button
                  type="button"
                  onClick={() => setSimilarityScore(Math.max(75, similarityScore - 1))}
                  aria-label="Decrease minimum code similarity"
                  className="w-8 h-8 rounded-md bg-slate-100 dark:bg-slate-700 hover:bg-slate-200 dark:hover:bg-slate-600 flex items-center justify-center cursor-pointer"
                >
                  −
                </button>
                <input
                  id="minimum-code-similarity"
                  type="range"
                  min="75"
                  max="100"
                  step="1"
                  value={similarityScore}
                  onChange={(e) => setSimilarityScore(Number(e.target.value))}
                  aria-label="Minimum code similarity"
                  aria-valuetext={`${similarityScore} percent`}
                  className="flex-1 cursor-pointer"
                />
                <button
                  type="button"
                  onClick={() => setSimilarityScore(Math.min(100, similarityScore + 1))}
                  aria-label="Increase minimum code similarity"
                  className="w-8 h-8 rounded-md bg-slate-100 dark:bg-slate-700 hover:bg-slate-200 dark:hover:bg-slate-600 flex items-center justify-center cursor-pointer"
                >
                  +
                </button>
              </div>
            </div>

            {/* Omnimodels */}
            <div>
              <h3 className="text-sm font-semibold text-slate-700 dark:text-slate-300 mb-3">
                Omnimodels
              </h3>
              <div className="space-y-2">
                <label className="flex items-center space-x-2 cursor-pointer">
                  <input
                    type="checkbox"
                    checked={omnimodel1}
                    onChange={(e) => setOmnimodel1(e.target.checked)}
                    className="rounded cursor-pointer"
                  />
                  <span className="text-sm text-slate-600 dark:text-slate-400">Omnimodel 1</span>
                </label>
                <label className="flex items-center space-x-2 cursor-pointer">
                  <input
                    type="checkbox"
                    checked={omnimodel2}
                    onChange={(e) => setOmnimodel2(e.target.checked)}
                    className="rounded cursor-pointer"
                  />
                  <span className="text-sm text-slate-600 dark:text-slate-400">Omnimodel 2</span>
                </label>
              </div>
            </div>

            {/* ICG Toggle */}
            <div>
              <div className="flex items-center justify-between">
                <span id="show-icg-entries-label" className="text-sm font-semibold text-slate-700 dark:text-slate-300">
                  Show ICG Entries
                </span>
                <button
                  type="button"
                  role="switch"
                  aria-checked={showICG}
                  aria-labelledby="show-icg-entries-label"
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
              </div>
            </div>

            {/* Selected Node Details */}
            {selectedNode && (
              <div className="border-t border-slate-200 dark:border-slate-700 pt-4">
                <div className="mb-3 flex items-center justify-between gap-2">
                  <h3 className="text-sm font-semibold text-slate-700 dark:text-slate-300">
                    Selected Node
                  </h3>
                  <button
                    type="button"
                    onClick={toggleInfoFreeze}
                    className={`text-xs font-medium px-2 py-1 rounded-md border transition-colors ${
                      isInfoFrozen
                        ? 'border-amber-400 text-amber-600 bg-amber-50 hover:bg-amber-100'
                        : 'border-slate-300 text-slate-600 hover:bg-slate-100'
                    }`}
                    title={isInfoFrozen ? 'Unlock info panel (L)' : 'Lock info panel (L)'}
                  >
                    {isInfoFrozen ? 'Unlock' : 'Lock'}
                  </button>
                </div>
                <div className="bg-slate-50 dark:bg-slate-900 rounded-lg p-3 text-sm">
                  {isInfoFrozen && (
                    <p className="mb-2 text-xs font-semibold uppercase tracking-wide text-amber-500">
                      Info locked
                    </p>
                  )}
                  <p className="text-slate-600 dark:text-slate-400">
                    <span className="font-medium">Name:</span> {selectedNodeDisplayName}
                  </p>
                  <p className="text-slate-600 dark:text-slate-400">
                    <span className="font-medium">Class:</span>{' '}
                    {formatIonClassForDisplay(selectedNode.original_model?.ion_class || selectedNode.ion_class)}
                  </p>
                  <p className="text-slate-600 dark:text-slate-400">
                    <span className="font-medium">ICG:</span> {selectedNode.icg ? 'Yes' : 'No'}
                  </p>
                  {(selectedNodeOmnimodel1 || selectedNodeOmnimodel2) && (
                    <p className="text-slate-600 dark:text-slate-400">
                      <span className="font-medium">Omnimodel:</span>{' '}
                      {[
                        selectedNodeOmnimodel1 ? '1' : null,
                        selectedNodeOmnimodel2 ? '2' : null,
                      ].filter(Boolean).join(' & ')}
                    </p>
                  )}
                  <button
                    type="button"
                    onClick={() => openSourceForNode(selectedNode)}
                    className="mt-3 inline-flex w-full items-center justify-center gap-2 rounded-md bg-slate-800 px-3 py-2 text-sm font-medium text-white transition-colors hover:bg-slate-700 focus-visible:outline-none focus-visible:ring-2 focus-visible:ring-blue-500 focus-visible:ring-offset-2 dark:bg-slate-200 dark:text-slate-900 dark:hover:bg-white"
                  >
                    View source code
                    <kbd aria-hidden="true" className="rounded bg-white/15 px-1.5 py-0.5 text-xs dark:bg-slate-900/10">S</kbd>
                  </button>
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
                  <span className="text-slate-600 dark:text-slate-400">Original/Ancestor Node</span>
                </div>
                <div className="pt-2 text-xs text-slate-500 dark:text-slate-500">
                  Node size represents number of identical models
                </div>
              </div>
            </div>

            {/* Keyboard Shortcuts */}
            <div className="border-t border-slate-200 dark:border-slate-700 pt-4">
              <h3 className="text-sm font-semibold text-slate-700 dark:text-slate-300 mb-3">
                Shortcuts
              </h3>
              <div className="space-y-1.5 text-xs text-slate-600 dark:text-slate-400">
                <div><kbd className="px-1.5 py-0.5 bg-slate-200 dark:bg-slate-700 rounded text-xs">Click</kbd> Select node</div>
                <div><kbd className="px-1.5 py-0.5 bg-slate-200 dark:bg-slate-700 rounded text-xs">Shift</kbd> + <kbd className="px-1.5 py-0.5 bg-slate-200 dark:bg-slate-700 rounded text-xs">Click</kbd> Add/remove from selection</div>
                <div><kbd className="px-1.5 py-0.5 bg-slate-200 dark:bg-slate-700 rounded text-xs">Drag</kbd> Box select</div>
                <div><kbd className="px-1.5 py-0.5 bg-slate-200 dark:bg-slate-700 rounded text-xs">Right-Click</kbd> Context menu</div>
                <div><kbd className="px-1.5 py-0.5 bg-slate-200 dark:bg-slate-700 rounded text-xs">A</kbd> Trace ancestors &nbsp; <kbd className="px-1.5 py-0.5 bg-slate-200 dark:bg-slate-700 rounded text-xs">S</kbd> Fetch source</div>
                <div><kbd className="px-1.5 py-0.5 bg-slate-200 dark:bg-slate-700 rounded text-xs">D</kbd> Diff selected &nbsp; <kbd className="px-1.5 py-0.5 bg-slate-200 dark:bg-slate-700 rounded text-xs">←</kbd><kbd className="px-1.5 py-0.5 bg-slate-200 dark:bg-slate-700 rounded text-xs">→</kbd> Step through diffs</div>
                <div><kbd className="px-1.5 py-0.5 bg-slate-200 dark:bg-slate-700 rounded text-xs">F</kbd> Fit to view &nbsp; <kbd className="px-1.5 py-0.5 bg-slate-200 dark:bg-slate-700 rounded text-xs">L</kbd> Lock panel &nbsp; <kbd className="px-1.5 py-0.5 bg-slate-200 dark:bg-slate-700 rounded text-xs">Esc</kbd> Close/clear</div>
              </div>
            </div>

            {/* Font Size Scaler */}
            <div className="border-t border-slate-200 dark:border-slate-700 pt-4">
              <div className="flex items-center justify-between">
                <span className="text-sm font-semibold text-slate-700 dark:text-slate-300">Font size</span>
                <div className="flex items-center gap-2">
                  <button
                    type="button"
                    onClick={() => {
                      setLeftFontScale((cur) => Math.max(0.8, Math.round((cur - 0.05) * 100) / 100));
                    }}
                    className="w-7 h-7 rounded-md border border-slate-300 dark:border-slate-600 text-slate-700 dark:text-slate-300 hover:bg-slate-100 dark:hover:bg-slate-700 flex items-center justify-center cursor-pointer"
                    aria-label="Decrease font size"
                  >
                    −
                  </button>
                  <span className="text-xs tabular-nums text-slate-600 dark:text-slate-400 select-none">
                    {`${Math.round(leftFontScale * 100)}%`}
                  </span>
                  <button
                    type="button"
                    onClick={() => {
                      setLeftFontScale((cur) => Math.min(1.5, Math.round((cur + 0.05) * 100) / 100));
                    }}
                    className="w-7 h-7 rounded-md border border-slate-300 dark:border-slate-600 text-slate-700 dark:text-slate-300 hover:bg-slate-100 dark:hover:bg-slate-700 flex items-center justify-center cursor-pointer"
                    aria-label="Increase font size"
                  >
                    +
                  </button>
                </div>
              </div>
            </div>
          </div>
        </div>

        {/* Visualization Area */}
        <div className="flex-1 flex flex-col min-w-0 relative">
          <div className="flex-1 bg-white dark:bg-slate-800 m-4 rounded-lg shadow-lg p-4 min-h-0 relative">
            {/* Info Freeze Indicator */}
            {isInfoFrozen && (
              <div className="absolute top-3 left-3 z-20 inline-flex items-center gap-2 rounded-full bg-amber-50 border border-amber-200 px-3 py-1 text-xs font-semibold text-amber-700 shadow-sm">
                <svg className="w-3.5 h-3.5" viewBox="0 0 24 24" fill="none" stroke="currentColor" strokeWidth="2">
                  <path strokeLinecap="round" strokeLinejoin="round" d="M9 3v2M15 3v2M7 7h10M5 11h14M6 21l-2-2m14 2 2-2M5 15h14M9 21v-2m6 2v-2" />
                </svg>
                Info locked (press L to unlock)
              </div>
            )}

            {isGroupedLayoutLoading && (
              <div role="status" className="absolute right-3 top-3 z-20 rounded-full border border-slate-200 bg-white/95 px-3 py-1 text-xs font-medium text-slate-600 shadow-sm dark:border-slate-700 dark:bg-slate-900/95 dark:text-slate-300">
                Loading grouped layout…
              </div>
            )}
            {groupedLayoutError && (
              <div role="alert" className="absolute left-3 right-3 top-3 z-20 rounded border border-red-200 bg-white p-3 text-sm text-red-700 dark:border-red-900 dark:bg-slate-900 dark:text-red-400">
                {groupedLayoutError}
                <button type="button" onClick={() => setLayoutRetry((value) => value + 1)} className="ml-2 underline">
                  Retry layout
                </button>
              </div>
            )}

            <svg 
              ref={svgRef} 
              className="w-full h-full border-2 border-slate-200 dark:border-slate-600 rounded"
              style={{ 
                visibility: groupedLayoutCount > 1 && !fixedLocationCircles[groupedLayoutCount] ? 'hidden' : 'visible',
                background: resolvedTheme === 'dark'
                  ? '#1e293b'
                  : '#fafafa'
              }}
            ></svg>
          </div>

          {/* Group summary overlay moved to right panel to avoid overlap */}
        </div>
        
        {/* Right Information Panel */}
        <div className="w-80 bg-white dark:bg-slate-800 border-l border-slate-200 dark:border-slate-700 overflow-y-auto scrollbar-none flex-shrink-0">
          <div className="p-4 space-y-4">
            {/* Groups Summary */}
            {groupSummaries.length > 0 && (
              <div className="bg-slate-50 dark:bg-slate-900 rounded-lg p-4">
                <h3 className="text-sm font-semibold text-amber-500 mb-3">Groups</h3>
                <div className="space-y-2 text-sm">
                  {groupSummaries.map((g, idx) => (
                    <div key={idx} className="flex items-baseline justify-between">
                      <div className="pr-2 text-slate-800 dark:text-slate-200">{g.key}</div>
                      <div className="text-slate-500 dark:text-slate-400 whitespace-nowrap">{g.nodeCount} nodes</div>
                    </div>
                  ))}
                </div>
              </div>
            )}

            {/* Selection Summary */}
            <div className="bg-slate-50 dark:bg-slate-900 rounded-lg p-4">
              <div className="flex items-center justify-between">
                <h3 className="text-sm font-semibold text-slate-700 dark:text-slate-300">Selection</h3>
                <div className="flex items-center gap-2">
                  <button
                    type="button"
                    onClick={() => {
                      const all = selectionBuckets.all;
                      if (all.length === 0) return;
                      const lines: string[] = [];
                      all.forEach((node) => {
                        const uniq = (node.original_model?.unique_modelDB_mod_id as string | undefined) || String(node.id);
                        lines.push(uniq);
                        const identicals = Array.isArray(node.identical_models) ? node.identical_models : [];
                        identicals.forEach((e) => {
                          const id = (e as any)?.unique_modelDB_mod_id as string | undefined;
                          if (id) lines.push(id);
                        });
                      });
                      const seen = new Set<string>();
                      const output = lines.filter((l) => (l && !seen.has(l) && seen.add(l)));
                      const text = output.join('\n');
                      navigator.clipboard.writeText(text).catch(() => {});
                    }}
                    className="px-2 py-1 text-xs rounded border border-slate-300 dark:border-slate-600 hover:bg-slate-100 dark:hover:bg-slate-700"
                  >Copy IDs</button>
                  <button
                    type="button"
                    onClick={() => {
                      const rows = selectionSummary.rows || [];
                      const esc = (v: unknown) => {
                        const s = (v ?? '').toString();
                        return /[",\n]/.test(s) ? '"' + s.replace(/"/g, '""') + '"' : s;
                      };
                      const csv: string[] = [];
                      csv.push(['node_id','unique_id','name','ion_class','copies'].map(esc).join(','));
                      rows.forEach((r) => {
                        const publicIonClass = r.ion_class === '—' ? '' : (r.ion_class || '');
                        csv.push([r.node_id, r.unique_id, r.name || '', publicIonClass, String(r.copies)].map(esc).join(','));
                      });
                      const blob = new Blob([csv.join('\n')], { type: 'text/csv;charset=utf-8' });
                      const url = URL.createObjectURL(blob);
                      const a = document.createElement('a');
                      a.href = url;
                      a.download = 'selection_summary.csv';
                      document.body.appendChild(a);
                      a.click();
                      document.body.removeChild(a);
                      URL.revokeObjectURL(url);
                    }}
                    className="px-2 py-1 text-xs rounded border border-slate-300 dark:border-slate-600 hover:bg-slate-100 dark:hover:bg-slate-700"
                  >Export CSV</button>
                </div>
              </div>
              <div className="mt-3 grid grid-cols-2 gap-3 text-sm">
                <div className="bg-white dark:bg-slate-800 rounded p-3 border border-slate-200 dark:border-slate-700">
                  <div className="text-xs uppercase tracking-wide text-slate-500">Selected nodes</div>
                  <div className="mt-1 text-lg font-semibold tabular-nums text-slate-900 dark:text-white">{selectionSummary.nodes}</div>
                </div>
                <div className="bg-white dark:bg-slate-800 rounded p-3 border border-slate-200 dark:border-slate-700">
                  <div className="text-xs uppercase tracking-wide text-slate-500">Total copies</div>
                  <div className="mt-1 text-lg font-semibold tabular-nums text-slate-900 dark:text-white">{selectionSummary.copies}</div>
                </div>
              </div>
              {selectionSummary.rows.length > 0 && (
                <div className="mt-3 max-h-48 overflow-y-auto overflow-x-hidden border border-slate-200 dark:border-slate-700 rounded">
                  <table className="w-full text-xs table-fixed">
                    <thead className="bg-slate-100 dark:bg-slate-800 text-slate-600 dark:text-slate-300">
                      <tr>
                        <th className="text-left px-2 py-1 w-[45%]">Node</th>
                        <th className="text-left px-2 py-1 w-[45%]">Unique ID</th>
                        <th className="text-right px-2 py-1 w-[10%]">Copies</th>
                      </tr>
                    </thead>
                    <tbody>
                      {selectionSummary.rows.map((r, i) => (
                        <tr key={`${r.node_id}-${i}`} className="odd:bg-white even:bg-slate-50 dark:odd:bg-slate-900 dark:even:bg-slate-800">
                          <td className="px-2 py-1 text-slate-700 dark:text-slate-300 truncate" title={r.name || r.node_id}>{r.name || r.node_id}</td>
                          <td className="px-2 py-1 text-slate-500 dark:text-slate-400 font-mono truncate" title={r.unique_id}>{r.unique_id}</td>
                          <td className="px-2 py-1 text-right text-slate-700 dark:text-slate-300 tabular-nums">{r.copies}</td>
                        </tr>
                      ))}
                    </tbody>
                  </table>
                </div>
              )}
            </div>
            {/* Node Details Box */}
            <div className="bg-slate-50 dark:bg-slate-900 rounded-lg p-4">
              <h3 className="text-sm font-semibold text-amber-500 mb-3">
                Node Details
              </h3>
              {selectedNode ? (
                <div className="space-y-2 text-sm">
                  <div>
                    <span className="text-amber-500 font-medium">ModelDB ID: </span>
                    <span className="text-slate-700 dark:text-slate-300">
                      {selectedNode.original_model?.unique_modelDB_mod_id || selectedNode.id}
                    </span>
                  </div>
                  <div>
                    <span className="text-amber-500 font-medium">Ion Class: </span>
                    <span className="text-slate-700 dark:text-slate-300">
                      {formatIonClassForDisplay(selectedNode.original_model?.ion_class || selectedNode.ion_class)}
                    </span>
                  </div>
                  <div>
                    <span className="text-amber-500 font-medium">Copies: </span>
                    <span className="text-slate-700 dark:text-slate-300">
                      {selectedNode.num_of_identicals || 1}
                    </span>
                  </div>
                  <div>
                    <span className="text-amber-500 font-medium">Year: </span>
                    <span className="text-slate-700 dark:text-slate-300">
                      {selectedNode.original_model?.Year || 'N/A'}
                    </span>
                  </div>
                  {selectedNodeAuthors && (
                    <div>
                      <span className="text-amber-500 font-medium">Authors: </span>
                      <span className="text-slate-700 dark:text-slate-300 text-xs">
                        {selectedNodeAuthors}
                      </span>
                    </div>
                  )}
                  {selectedNode.original_model?.ICG && (
                    <div>
                      <span className="text-amber-500 font-medium">ICG Entry: </span>
                      <span className="text-green-500">Yes</span>
                    </div>
                  )}
                  {(selectedNodeOmnimodel1 || selectedNodeOmnimodel2) && (
                    <div>
                      <span className="text-amber-500 font-medium">Omnimodel: </span>
                      <span className="text-slate-700 dark:text-slate-300">
                        {[
                          selectedNodeOmnimodel1 ? '1' : null,
                          selectedNodeOmnimodel2 ? '2' : null,
                        ].filter(Boolean).join(' & ')}
                      </span>
                    </div>
                  )}
                  {selectedNodeOmnimodelHref && (selectedNodeOmnimodel1 || selectedNodeOmnimodel2) && (
                    <div>
                      <a
                        href={selectedNodeOmnimodelHref}
                        target="_blank"
                        rel="noopener noreferrer"
                        className="inline-flex items-center gap-1 text-sm font-medium text-blue-600 hover:text-blue-700"
                      >
                        View omnimodel markdown
                        <svg className="w-4 h-4" fill="none" stroke="currentColor" viewBox="0 0 24 24">
                          <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2} d="M17 7l-10 10M8 7h9v9" />
                        </svg>
                      </a>
                    </div>
                  )}
                </div>
              ) : (
                <p className="text-sm text-slate-500 dark:text-slate-400 italic">
                  Click on a node to view details
                </p>
              )}
            </div>

            {/* Copies Box */}
            <div className="bg-slate-50 dark:bg-slate-900 rounded-lg p-4">
              <h3 className="text-sm font-semibold text-amber-500 mb-3">Copies</h3>
              {selectedNode?.identical_models && selectedNode.identical_models.length > 0 ? (
                <div className="max-h-48 overflow-y-auto scrollbar-none space-y-1">
                  {selectedNode.identical_models.map((model: any, idx: number) => (
                    <div key={idx} className="text-sm text-slate-700 dark:text-slate-300">
                      {model.unique_modelDB_mod_id || model}
                    </div>
                  ))}
                </div>
              ) : (
                <p className="text-sm text-slate-500 dark:text-slate-400 italic">
                  {selectedNode ? 'No identical models' : 'Select a node to view identical models'}
                </p>
              )}
            </div>

            {/* Removed right-panel prints: statistics and selection info */}
          </div>
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
              {getNodeDisplayName(contextMenu.node)}
            </div>
          <div className="text-xs text-slate-500 dark:text-slate-400">
              {formatIonClassForDisplay(contextMenu.node.original_model?.ion_class || contextMenu.node.ion_class)} • {contextMenu.node.icg ? 'ICG' : 'ModelDB'}
          </div>
          </div>
          
          {contextMenuOmnimodelHref && (contextMenuOmnimodel1 || contextMenuOmnimodel2) && (
            <a
              href={contextMenuOmnimodelHref}
              target="_blank"
              rel="noopener noreferrer"
              onClick={() => setContextMenu(null)}
              className="w-full px-4 py-2 text-left text-sm text-slate-700 dark:text-slate-300 hover:bg-slate-100 dark:hover:bg-slate-700 flex items-center space-x-2"
            >
              <svg className="w-4 h-4" fill="none" stroke="currentColor" viewBox="0 0 24 24">
                <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2} d="M17 7l-10 10M8 7h9v9" />
              </svg>
              <span>Open omnimodel markdown</span>
            </a>
          )}

          <button
            type="button"
            className="w-full px-4 py-2 text-left text-sm text-slate-700 dark:text-slate-300 hover:bg-slate-100 dark:hover:bg-slate-700 flex items-center space-x-2 cursor-pointer"
            onClick={() => openSourceForNode(contextMenu.node)}
          >
            <svg className="w-4 h-4" fill="none" stroke="currentColor" viewBox="0 0 24 24" aria-hidden="true">
              <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2} d="M8 9l3 3-3 3m5 0h3M6 3h12a2 2 0 012 2v14a2 2 0 01-2 2H6a2 2 0 01-2-2V5a2 2 0 012-2z" />
            </svg>
            <span>View source code</span>
          </button>

          <button
            type="button"
            className="w-full px-4 py-2 text-left text-sm text-slate-700 dark:text-slate-300 hover:bg-slate-100 dark:hover:bg-slate-700 flex items-center space-x-2 cursor-pointer"
            onClick={async () => {
              const combined = new Map<string, Node>();
              selectionBuckets.all.forEach((node) => {
                combined.set(String(node.id), node);
              });
              combined.set(String(contextMenu.node.id), contextMenu.node);
              await generateAwesomeDiffs(Array.from(combined.values()), false);
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
              const nodeId = String(contextMenu.node.id);
              commitSelectionState((prev) => {
                const selectedSet = new Set(prev.selectedIds);
                const sourceSet = new Set(prev.sourceIds);
                let anchorId = prev.anchorId;

                if (selectedSet.has(nodeId)) {
                  selectedSet.delete(nodeId);
                  if (anchorId === nodeId) {
                    anchorId = null;
                  }
                } else if (sourceSet.has(nodeId)) {
                  sourceSet.delete(nodeId);
                  selectedSet.add(nodeId);
                  anchorId = nodeId;
                } else {
                  selectedSet.add(nodeId);
                  anchorId = nodeId;
                }

                return { selectedIds: Array.from(selectedSet), sourceIds: Array.from(sourceSet), anchorId };
              });
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
              const nodeId = String(contextMenu.node.id);
              commitSelectionState((prev) => {
                const selectedSet = new Set(prev.selectedIds);
                const sourceSet = new Set(prev.sourceIds);
                let anchorId = prev.anchorId;

                if (sourceSet.has(nodeId)) {
                  sourceSet.delete(nodeId);
                } else {
                  sourceSet.add(nodeId);
                  if (selectedSet.has(nodeId)) {
                    selectedSet.delete(nodeId);
                  }
                  if (!anchorId) {
                    anchorId = nodeId;
                  }
                }

                return { selectedIds: Array.from(selectedSet), sourceIds: Array.from(sourceSet), anchorId };
              });
            }}
          >
            <svg className="w-4 h-4" fill="none" stroke="currentColor" viewBox="0 0 24 24">
              <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2} d="M11.049 2.927c.3-.921 1.603-.921 1.902 0l1.519 4.674a1 1 0 00.95.69h4.915c.969 0 1.371 1.24.588 1.81l-3.976 2.888a1 1 0 00-.363 1.118l1.518 4.674c.3.922-.755 1.688-1.538 1.118l-3.976-2.888a1 1 0 00-1.176 0l-3.976 2.888c-.783.57-1.838-.197-1.538-1.118l1.518-4.674a1 1 0 00-.363-1.118l-3.976-2.888c-.784-.57-.38-1.81.588-1.81h4.914a1 1 0 00.951-.69l1.519-4.674z" />
            </svg>
            <span>Toggle as source</span>
          </button>
        </div>
      )}

      {sourceView && (
        <div className="fixed inset-0 z-[60] flex items-center justify-center bg-black/75 p-3 sm:p-6">
          <div
            ref={sourceDialogRef}
            role="dialog"
            aria-modal="true"
            aria-labelledby="source-code-dialog-title"
            aria-describedby="source-code-dialog-description"
            className="flex h-[92vh] w-full max-w-6xl flex-col overflow-hidden rounded-xl border border-slate-200 bg-white shadow-2xl dark:border-slate-700 dark:bg-slate-900"
            onKeyDown={(event) => {
              if (event.key !== 'Tab') return;
              const focusable = sourceDialogRef.current?.querySelectorAll<HTMLElement>(
                'button:not([disabled]), [href], input:not([disabled]), select:not([disabled]), textarea:not([disabled]), [tabindex]:not([tabindex="-1"])'
              );
              if (!focusable || focusable.length === 0) return;
              const first = focusable[0];
              const last = focusable[focusable.length - 1];
              if (event.shiftKey && document.activeElement === first) {
                event.preventDefault();
                last.focus();
              } else if (!event.shiftKey && document.activeElement === last) {
                event.preventDefault();
                first.focus();
              }
            }}
          >
            <div className="flex items-start justify-between gap-4 border-b border-slate-700 bg-slate-800 px-5 py-4 text-white sm:px-6">
              <div className="min-w-0">
                <h2 id="source-code-dialog-title" className="text-xl font-semibold sm:text-2xl">
                  Source code
                </h2>
                <p id="source-code-dialog-description" className="mt-1 truncate font-mono text-sm text-slate-300">
                  {getNodeUniqueId(sourceView.node) || getNodeDisplayName(sourceView.node)}
                </p>
              </div>
              <button
                ref={sourceDialogCloseRef}
                type="button"
                onClick={closeSourceView}
                aria-label="Close source code viewer"
                className="rounded-lg p-2 text-slate-200 transition-colors hover:bg-white/10 hover:text-white focus-visible:outline-none focus-visible:ring-2 focus-visible:ring-white"
              >
                <svg className="h-6 w-6" fill="none" stroke="currentColor" viewBox="0 0 24 24" aria-hidden="true">
                  <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2} d="M6 18L18 6M6 6l12 12" />
                </svg>
              </button>
            </div>

            <div className="min-h-0 flex-1 bg-slate-950 p-3 sm:p-5">
              {sourceView.status === 'loading' && (
                <div role="status" aria-live="polite" className="flex h-full items-center justify-center text-slate-300">
                  <div className="text-center">
                    <div className="mx-auto mb-4 h-10 w-10 animate-spin rounded-full border-2 border-slate-600 border-t-blue-400" />
                    <p>Loading source code…</p>
                  </div>
                </div>
              )}

              {sourceView.status === 'error' && (
                <div role="alert" className="flex h-full items-center justify-center">
                  <div className="max-w-lg rounded-lg border border-red-800 bg-red-950/60 p-6 text-center text-red-100">
                    <h3 className="text-lg font-semibold">Source code unavailable</h3>
                    <p className="mt-2 text-sm text-red-200">{sourceView.error}</p>
                    <button
                      type="button"
                      onClick={() => openSourceForNode(sourceView.node)}
                      className="mt-5 rounded-md bg-red-100 px-4 py-2 text-sm font-semibold text-red-950 transition-colors hover:bg-white focus-visible:outline-none focus-visible:ring-2 focus-visible:ring-white"
                    >
                      Try again
                    </button>
                  </div>
                </div>
              )}

              {sourceView.status === 'ready' && (
                <pre
                  tabIndex={0}
                  aria-label={`Source code for ${getNodeUniqueId(sourceView.node) || getNodeDisplayName(sourceView.node)}`}
                  className="h-full overflow-auto rounded-lg border border-slate-800 bg-slate-950 p-4 font-mono text-[13px] leading-6 text-slate-100 outline-none focus-visible:ring-2 focus-visible:ring-blue-500"
                ><code>{sourceView.sourceCode}</code></pre>
              )}
            </div>

            <div className="flex items-center justify-between gap-4 border-t border-slate-200 bg-slate-50 px-5 py-3 text-xs text-slate-600 dark:border-slate-700 dark:bg-slate-800 dark:text-slate-300 sm:px-6">
              <span>Source is loaded read-only from the selected model.</span>
              <span className="whitespace-nowrap"><kbd className="rounded bg-slate-200 px-1.5 py-0.5 dark:bg-slate-700">Esc</kbd> Close</span>
            </div>
          </div>
        </div>
      )}

      {/* Awesome Diff View Overlay */}
      {showDiffView && diffCombinations.length > 0 && (
        <div className="fixed inset-0 bg-black bg-opacity-75 z-50 flex items-center justify-center p-2">
          <div className="bg-white dark:bg-slate-800 rounded-xl shadow-2xl w-[98%] max-w-[1920px] h-[98vh] overflow-hidden border border-slate-200 dark:border-slate-700 flex flex-col">
            {/* Awesome Diff Header */}
            <div className="bg-slate-800 text-white p-6">
              <div className="flex items-center justify-between">
                <div className="flex items-center space-x-4">
                  <div className="bg-white bg-opacity-20 rounded-lg p-3">
                    <span className="text-2xl">📄</span>
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
                    inFlightDiffKeysRef.current.clear();
                  }}
                  className="bg-white bg-opacity-20 hover:bg-opacity-30 rounded-lg p-2 transition-all duration-200 cursor-pointer"
                >
                  <svg className="w-6 h-6 text-white" fill="none" stroke="currentColor" viewBox="0 0 24 24" xmlns="http://www.w3.org/2000/svg">
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
                    <div className="w-3 h-3 rounded-full bg-amber-500"></div>
                    <span className="text-sm font-medium text-slate-700 dark:text-slate-300">Source:</span>
                    <span className="text-sm text-slate-600 dark:text-slate-400 font-mono">
                      {diffCombinations[currentCombinationIndex]?.source.name || diffCombinations[currentCombinationIndex]?.source.id}
                    </span>
                    {currentSourceIonClass && currentSourceIonClass !== '—' && (
                      <span className="text-xs text-slate-500 dark:text-slate-500 bg-slate-200 dark:bg-slate-700 px-2 py-1 rounded">
                        {currentSourceIonClass}
                      </span>
                    )}
                  </div>
                  <div className="text-slate-400">→</div>
                  <div className="flex items-center space-x-2">
                    <div className="w-3 h-3 rounded-full bg-blue-500"></div>
                    <span className="text-sm font-medium text-slate-700 dark:text-slate-300">Target:</span>
                    <span className="text-sm text-slate-600 dark:text-slate-400 font-mono">
                      {diffCombinations[currentCombinationIndex]?.target.name || diffCombinations[currentCombinationIndex]?.target.id}
                    </span>
                    {currentTargetIonClass && currentTargetIonClass !== '—' && (
                      <span className="text-xs text-slate-500 dark:text-slate-500 bg-slate-200 dark:bg-slate-700 px-2 py-1 rounded">
                        {currentTargetIonClass}
                      </span>
                    )}
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
                      ◀
                    </button>
                    <button
                      onClick={() => setCurrentCombinationIndex(Math.min(diffCombinations.length - 1, currentCombinationIndex + 1))}
                      disabled={currentCombinationIndex === diffCombinations.length - 1}
                      className="p-2 rounded-lg bg-slate-200 dark:bg-slate-700 hover:bg-slate-300 dark:hover:bg-slate-600 disabled:opacity-50 disabled:cursor-not-allowed transition-all duration-200"
                    >
                      ▶
                    </button>
                  </div>
                </div>
              </div>
            </div>

            {/* Progress Bar */}
            <div className="bg-slate-100 dark:bg-slate-800 h-1">
              <div 
                className="h-full bg-blue-500 transition-all duration-300 ease-out"
                style={{ width: `${((currentCombinationIndex + 1) / diffCombinations.length) * 100}%` }}
              />
            </div>

            {/* Diff Content */}
            <div className="flex-1 overflow-hidden">
              {diffCombinations[currentCombinationIndex]?.html ? (
                <div className="h-full p-4">
                  <iframe
                    title="Code diff"
                    className="w-full h-full rounded-md border border-slate-200 dark:border-slate-700 bg-white"
                    sandbox=""
                    referrerPolicy="no-referrer"
                    // Render escaped diff markup in a unique-origin sandbox.
                    srcDoc={diffCombinations[currentCombinationIndex].html}
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

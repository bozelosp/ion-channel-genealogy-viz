import fs from 'fs';
import path from 'path';

import { SHOW_OTHER_ION_CLASS, shouldSuppressIonClass } from './ionClasses';

interface RawNode {
  id: string | number;
  ion_class?: string | null;
  icg?: boolean;
  omnimodel1?: unknown;
  omnimodel2?: unknown;
  original_model?: {
    ion_class?: string | null;
    ICG?: boolean;
    unique_modelDB_mod_id?: string | null;
    [key: string]: unknown;
  };
  identical_models?: Array<{
    unique_modelDB_mod_id?: string | null;
    [key: string]: unknown;
  }>;
}

interface RawLink {
  source: string | number | { id: string | number };
  target: string | number | { id: string | number };
}

interface RawNetworkData {
  nodes: RawNode[];
  links: RawLink[];
}

interface AggregatedCounts {
  nodes: number;
  links: number;
  uniqueIds: number;
  gatingUniqueIds: number;
  fittedUniqueIds: number;
}

interface NetworkAggregates {
  total: AggregatedCounts;
  filtered: AggregatedCounts;
}

interface Metrics {
  filesProcessed: number;
  uniqueImplementations: number;
  similarityEdges: number;
  percentDuplicates: number;
  channelsWithGates: number;
  channelsWithGatesPercent: number;
  omnimodels: number;
  omnimodelCoveragePercent: number;
}

const NETWORK_DATA_PATH = path.resolve(process.cwd(), 'public', 'data', 'network_data.json');
const OMNIMODEL_ROOT = path.resolve(process.cwd(), 'public', 'data', 'omnimodels');

function readNetworkData(): RawNetworkData | null {
  try {
    const content = fs.readFileSync(NETWORK_DATA_PATH, 'utf-8');
    const parsed = JSON.parse(content) as RawNetworkData;
    if (!Array.isArray(parsed?.nodes) || !Array.isArray(parsed?.links)) {
      return null;
    }
    return parsed;
  } catch {
    return null;
  }
}

function toNodeId(value: string | number | { id?: unknown } | undefined): string | null {
  if (value === null || value === undefined) return null;
  if (typeof value === 'string' || typeof value === 'number') {
    return String(value);
  }
  if (typeof value === 'object' && 'id' in value) {
    const raw = (value as { id?: unknown }).id;
    if (raw === null || raw === undefined) return null;
    return String(raw);
  }
  return null;
}

function collectUniqueIds(node: RawNode): Set<string> {
  const ids = new Set<string>();
  const primary = node.original_model?.unique_modelDB_mod_id;
  if (typeof primary === 'string' && primary.length > 0) {
    ids.add(primary);
  }

  if (Array.isArray(node.identical_models)) {
    for (const ident of node.identical_models) {
      const candidate = ident?.unique_modelDB_mod_id;
      if (typeof candidate === 'string' && candidate.length > 0) {
        ids.add(candidate);
      }
    }
  }

  return ids;
}

function toBoolean(value: unknown): boolean {
  if (typeof value === 'boolean') return value;
  if (typeof value === 'number') return value !== 0;
  if (typeof value === 'string') {
    const lowered = value.trim().toLowerCase();
    return lowered === 'true' || lowered === '1' || lowered === 'yes';
  }
  return false;
}

function hasOmnimodelFlag(node: RawNode): boolean {
  const orig = node.original_model ?? {};
  return (
    toBoolean(node.omnimodel1) ||
    toBoolean(node.omnimodel2) ||
    toBoolean(orig['Supermodel 1']) ||
    toBoolean(orig['Supermodel 2'])
  );
}

function hasGatingFlag(node: RawNode): boolean {
  return Boolean(node.icg || node.original_model?.ICG);
}

function resolveIonClass(node: RawNode): string | null {
  return (
    (node.original_model?.ion_class ?? node.ion_class ?? null) as string | null
  );
}

function aggregateNetwork(data: RawNetworkData): NetworkAggregates {
  const uniqueAll = new Set<string>();
  const uniqueFiltered = new Set<string>();
  const gatingAll = new Set<string>();
  const gatingFiltered = new Set<string>();
  const fittedAll = new Set<string>();
  const fittedFiltered = new Set<string>();
  const filteredNodeIds = new Set<string>();

  let filteredNodeCount = 0;

  for (const node of data.nodes) {
    const ionClass = resolveIonClass(node);
    const suppressed = shouldSuppressIonClass(ionClass ?? undefined);
    const nodeIds = collectUniqueIds(node);
    const hasGating = hasGatingFlag(node);
    const hasOmni = hasOmnimodelFlag(node);

    for (const id of nodeIds) {
      uniqueAll.add(id);
      if (!suppressed) {
        uniqueFiltered.add(id);
      }
      if (hasGating) {
        gatingAll.add(id);
        if (!suppressed) gatingFiltered.add(id);
      }
      if (hasGating && hasOmni) {
        fittedAll.add(id);
        if (!suppressed) fittedFiltered.add(id);
      }
    }

    if (!suppressed) {
      filteredNodeCount += 1;
      const rawId = toNodeId(node.id);
      if (rawId) {
        filteredNodeIds.add(rawId);
      }
    }
  }

  let filteredLinkCount = 0;
  for (const link of data.links) {
    const sourceId = toNodeId(link.source);
    const targetId = toNodeId(link.target);
    if (!sourceId || !targetId) {
      continue;
    }
    if (filteredNodeIds.has(sourceId) && filteredNodeIds.has(targetId)) {
      filteredLinkCount += 1;
    }
  }

  const total: AggregatedCounts = {
    nodes: data.nodes.length,
    links: data.links.length,
    uniqueIds: uniqueAll.size,
    gatingUniqueIds: gatingAll.size,
    fittedUniqueIds: fittedAll.size,
  };

  const filtered: AggregatedCounts = {
    nodes: filteredNodeCount,
    links: filteredLinkCount,
    uniqueIds: uniqueFiltered.size,
    gatingUniqueIds: gatingFiltered.size,
    fittedUniqueIds: fittedFiltered.size,
  };

  return { total, filtered };
}

function safeReadDirectory(dirPath: string): string[] {
  try {
    return fs.readdirSync(dirPath);
  } catch {
    return [];
  }
}

function countOmnimodelDirectories(): { total: number; filtered: number } {
  const channels = safeReadDirectory(OMNIMODEL_ROOT);
  if (channels.length === 0) {
    return { total: 0, filtered: 0 };
  }

  let total = 0;
  let filtered = 0;

  for (const channel of channels) {
    const channelPath = path.join(OMNIMODEL_ROOT, channel);
    if (!fs.existsSync(channelPath) || !fs.statSync(channelPath).isDirectory()) {
      continue;
    }

    const modelDirs = safeReadDirectory(channelPath).filter((entry) => {
      const modelPath = path.join(channelPath, entry);
      return fs.existsSync(modelPath) && fs.statSync(modelPath).isDirectory();
    });

    total += modelDirs.length;
    if (!shouldSuppressIonClass(channel)) {
      filtered += modelDirs.length;
    }
  }

  return { total, filtered };
}

function computePercent(part: number, whole: number): number {
  if (!whole) return 0;
  return (part / whole) * 100;
}

function roundPercent(value: number): number {
  return Math.round(value * 10) / 10;
}

function buildMetrics(): Metrics {
  const data = readNetworkData();
  if (!data) {
    return {
      filesProcessed: 0,
      uniqueImplementations: 0,
      similarityEdges: 0,
      percentDuplicates: 0,
      channelsWithGates: 0,
      channelsWithGatesPercent: 0,
      omnimodels: 0,
      omnimodelCoveragePercent: 0,
    };
  }

  const aggregates = aggregateNetwork(data);
  const omniDirCounts = countOmnimodelDirectories();

  const scope = SHOW_OTHER_ION_CLASS ? aggregates.total : aggregates.filtered;
  const rawOmniCount = SHOW_OTHER_ION_CLASS ? omniDirCounts.total : omniDirCounts.filtered;
  const resolvedOmniCount = rawOmniCount > 0 ? Math.min(rawOmniCount, scope.fittedUniqueIds) : scope.fittedUniqueIds;

  const duplicateShare = scope.uniqueIds === 0
    ? 0
    : Math.max(0, ((scope.uniqueIds - scope.nodes) / scope.uniqueIds) * 100);

  const channelsWithGatesPercent = roundPercent(computePercent(scope.gatingUniqueIds, scope.uniqueIds));
  const omnimodelCoveragePercent = roundPercent(
    scope.gatingUniqueIds === 0 ? 0 : computePercent(resolvedOmniCount, scope.gatingUniqueIds),
  );

  return {
    filesProcessed: scope.uniqueIds,
    uniqueImplementations: scope.nodes,
    similarityEdges: scope.links,
    percentDuplicates: roundPercent(duplicateShare),
    channelsWithGates: scope.gatingUniqueIds,
    channelsWithGatesPercent,
    omnimodels: resolvedOmniCount,
    omnimodelCoveragePercent,
  };
}

export const metrics = buildMetrics();

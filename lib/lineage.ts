// Lineage of one .mod file family, reconstructed from the source-similarity
// graph. Each dated file is linked to its most similar earlier file in the
// same family; a file with no earlier relative starts its own lineage.
// Undated files and files with no linked relative are left out. Invalid or
// unexpected input throws, so a changed dataset fails the build instead of
// rendering a different lineage.

import { isModelSourceId } from './comparison';
import { graphIonClass, shouldSuppressIonClass } from './ionClasses';

export interface LineageNode {
  id: string;
  modelDb: string;
  name: string;
  year: number;
  /** Index of the most similar earlier file, or -1 for a lineage origin. */
  parent: number;
  /** Index of this node's lineage origin (itself for an origin). */
  origin: number;
  simParent: number | null;
  /** Similarity to the lineage origin, or null when they are not linked. */
  simOrigin: number | null;
  copies: number;
  descendants: number;
  hasChildren: boolean;
  visibleInGraph: boolean;
}

export interface Lineage {
  file: string;
  nodes: LineageNode[];
  /** Lowest link weight in the graph: links exist only at or above it. */
  linkThreshold: number;
}

export interface LineageLayout {
  size: number;
  points: { x: number; y: number }[];
  /** SVG path of the edge into each node from its parent ('' for origins). */
  edges: string[];
  rings: { year: number; r: number }[];
}

interface Member {
  key: number;
  id: string;
  modelDb: string;
  variant: number;
  year: number;
  copies: number;
  parent: number | null;
  simParent: number | null;
  children: number[];
  visibleInGraph: boolean;
}

function record(value: unknown): Record<string, unknown> | null {
  return value && typeof value === 'object' && !Array.isArray(value) ? (value as Record<string, unknown>) : null;
}

function byYearThenId(a: Member, b: Member): number {
  return a.year - b.year || a.id.localeCompare(b.id, 'en');
}

export function buildLineage(data: unknown, file: string, primaryOrigin: string): Lineage {
  const graph = record(data);
  if (!graph || !Array.isArray(graph.nodes) || !Array.isArray(graph.links)) {
    throw new Error('Network data must contain node and link arrays');
  }

  const family = new Map<number, Member>();
  for (const raw of graph.nodes) {
    const node = record(raw);
    const model = record(node?.original_model);
    if (!node || !model || model.mod_filename !== file) continue;
    const year = model.Year;
    if (typeof year !== 'number' || !Number.isInteger(year) || year < 1900) continue;
    const id = model.unique_modelDB_mod_id;
    const match = typeof id === 'string' && isModelSourceId(id) ? /^(\d+)_.+-ID-(\d+)$/.exec(id) : null;
    const copies = node.num_of_identicals;
    if (!match || typeof node.id !== 'number' || typeof copies !== 'number' || !Number.isInteger(copies) || copies < 0) {
      throw new Error(`Invalid ${file} node in network data`);
    }
    family.set(node.id, {
      key: node.id, id: id as string, modelDb: match[1], variant: Number(match[2]), year, copies,
      parent: null, simParent: null, children: [],
      visibleInGraph: !shouldSuppressIonClass(graphIonClass(model.ion_class, node.ion_class, node.ion_channel_class)),
    });
  }

  let linkThreshold = Infinity;
  const weights = new Map<number, Map<number, number>>();
  for (const raw of graph.links) {
    const link = record(raw);
    const weight = link?.weight;
    if (!link || typeof weight !== 'number' || !(weight > 0 && weight <= 100)) {
      throw new Error('Invalid link weight in network data');
    }
    linkThreshold = Math.min(linkThreshold, weight);
    const { source, target } = link;
    if (typeof source !== 'number' || typeof target !== 'number') continue;
    if (!family.has(source) || !family.has(target)) continue;
    for (const [a, b] of [[source, target], [target, source]]) {
      if (!weights.has(a)) weights.set(a, new Map());
      weights.get(a)!.set(b, weight);
    }
  }

  for (const member of family.values()) {
    let best: Member | null = null;
    let bestWeight = -Infinity;
    for (const [key, weight] of weights.get(member.key) ?? []) {
      const other = family.get(key)!;
      if (other.year >= member.year) continue;
      if (weight > bestWeight || (weight === bestWeight && best && byYearThenId(other, best) < 0)) {
        best = other;
        bestWeight = weight;
      }
    }
    if (best) {
      member.parent = best.key;
      member.simParent = bestWeight;
      best.children.push(member.key);
    }
  }

  const ordered: Member[] = [];
  const visit = (member: Member) => {
    ordered.push(member);
    member.children.map((key) => family.get(key)!).sort(byYearThenId).forEach(visit);
  };
  [...family.values()]
    .filter((member) => member.parent === null && member.children.length > 0)
    .sort(byYearThenId)
    .forEach(visit);

  if (ordered[0]?.id !== primaryOrigin) {
    throw new Error(`Expected ${primaryOrigin} to be the earliest ${file} lineage origin`);
  }

  const index = new Map(ordered.map((member, i) => [member.key, i]));
  const variantsPerEntry = new Map<string, number>();
  for (const member of ordered) variantsPerEntry.set(member.modelDb, (variantsPerEntry.get(member.modelDb) ?? 0) + 1);
  const nodes: LineageNode[] = [];
  ordered.forEach((member, i) => {
    const parent = member.parent === null ? -1 : index.get(member.parent)!;
    const origin = parent < 0 ? i : nodes[parent].origin;
    nodes.push({
      id: member.id,
      modelDb: member.modelDb,
      name: member.modelDb + (variantsPerEntry.get(member.modelDb)! > 1 ? '·' + member.variant : ''),
      year: member.year,
      parent,
      origin,
      simParent: member.simParent,
      simOrigin: origin === i ? null : weights.get(member.key)?.get(ordered[origin].key) ?? null,
      copies: member.copies,
      descendants: 0,
      hasChildren: member.children.length > 0,
      visibleInGraph: member.visibleInGraph,
    });
  });
  for (let i = nodes.length - 1; i >= 0; i -= 1) {
    const parent = nodes[i].parent;
    if (parent >= 0) nodes[parent].descendants += 1 + nodes[i].descendants;
  }

  return { file, nodes, linkThreshold };
}

const round = (value: number) => Math.round(value * 10) / 10;

/**
 * Radial layout: the primary origin sits at the centre, radius encodes the
 * year, and leaves are spread evenly by angle in depth-first order, with one
 * empty slot before each further lineage.
 */
export function layoutLineage({ nodes }: Lineage, size = 540): LineageLayout {
  const centre = size / 2;
  const firstYear = nodes[0].year;
  const lastYear = Math.max(...nodes.map((node) => node.year));
  const innerRadius = 22;
  const outerRadius = centre - 27;
  const radius = (year: number) =>
    innerRadius + ((year - firstYear) * (outerRadius - innerRadius)) / Math.max(1, lastYear - firstYear);

  const slots = new Array<number>(nodes.length).fill(0);
  let next = 0;
  nodes.forEach((node, i) => {
    if (i > 0 && node.parent < 0) next += 1;
    if (!node.hasChildren) slots[i] = next++;
  });
  const childSlots = nodes.map(() => [] as number[]);
  for (let i = nodes.length - 1; i >= 0; i -= 1) {
    if (nodes[i].hasChildren) {
      slots[i] = childSlots[i].reduce((sum, slot) => sum + slot, 0) / childSlots[i].length;
    }
    if (nodes[i].parent >= 0) childSlots[nodes[i].parent].push(slots[i]);
  }

  const step = (2 * Math.PI) / Math.max(1, next);
  const startAngle = -Math.PI / 2 - 0.25;
  const polar = nodes.map((node, i) =>
    i === 0 ? { r: 0, a: 0 } : { r: radius(node.year), a: startAngle + slots[i] * step },
  );
  const points = polar.map(({ r, a }) => ({ x: round(centre + r * Math.cos(a)), y: round(centre + r * Math.sin(a)) }));

  const edges = nodes.map((node, i) => {
    if (node.parent < 0) return '';
    const to = points[i];
    if (node.parent === 0) return `M${centre} ${centre}L${to.x} ${to.y}`;
    const from = points[node.parent];
    const { r, a } = polar[node.parent];
    const { a: b } = polar[i];
    let path = `M${from.x} ${from.y}`;
    if (Math.abs(b - a) > 1e-6) {
      path += `A${round(r)} ${round(r)} 0 0 ${b > a ? 1 : 0} ${round(centre + r * Math.cos(b))} ${round(centre + r * Math.sin(b))}`;
    }
    return `${path}L${to.x} ${to.y}`;
  });

  const rings: LineageLayout['rings'] = [];
  for (let year = Math.ceil((firstYear + 1) / 5) * 5; year <= lastYear; year += 5) {
    const r = radius(year);
    if (r >= 60) rings.push({ year, r: round(r) });
  }

  return { size, points, edges, rings };
}

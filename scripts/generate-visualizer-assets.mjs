#!/usr/bin/env node

import path from 'node:path';
import { promises as fs } from 'node:fs';
import { fileURLToPath } from 'node:url';
import { assertWritableWithin } from './path-safety.mjs';

const REPO_ROOT = path.resolve(path.dirname(fileURLToPath(import.meta.url)), '..');
const DATA_ROOT = path.join(REPO_ROOT, 'public', 'data');
const NETWORK_SOURCE = path.join(DATA_ROOT, 'network_data.json');
const INITIAL_NETWORK_OUTPUT = path.join(DATA_ROOT, 'network_data_initial.json');
const LAYOUT_SOURCE = path.join(DATA_ROOT, 'fixed_location_circles.json');
const LAYOUT_OUTPUT_ROOT = path.join(DATA_ROOT, 'fixed-location-circles');
const INITIAL_MINIMUM_COPIES = 2;
const MINIMUM_LAYOUT_GROUPS = 2;
const MAXIMUM_LAYOUT_GROUPS = 25;

function identifier(value) {
  if (typeof value === 'string' && value.length > 0) return `s:${value}`;
  if (typeof value === 'number' && Number.isFinite(value)) return `n:${value}`;
  throw new Error(`Invalid network identifier: ${JSON.stringify(value)}`);
}

const CHANNEL_ALIASES = new Map([
  ['na', 'Na'],
  ['k', 'K'],
  ['ca', 'Ca'],
  ['ih', 'IH'],
  ['h', 'IH'],
  ['kca', 'KCa'],
]);

function channelName(value) {
  if (typeof value !== 'string') return null;
  return CHANNEL_ALIASES.get(value.trim().toLowerCase()) ?? null;
}

function locationFromPath(value) {
  if (typeof value !== 'string') return null;
  const segments = value.split('/').filter(Boolean);
  const channelSegment = segments.find((segment) => segment.startsWith('icg-channels-'));
  if (!channelSegment) return null;
  const channel = channelName(channelSegment.slice('icg-channels-'.length));
  const slug = segments[segments.indexOf(channelSegment) + 1];
  return channel && slug ? [channel, slug] : null;
}

function omnimodelLocation(model) {
  if (Array.isArray(model?.ICG_entries)) {
    for (const entry of model.ICG_entries) {
      const fromInfo = locationFromPath(entry?.info?.mod_filepath ?? entry?.info?.modFilename);
      if (fromInfo) return fromInfo;
      const fromEntry = locationFromPath(entry?.mod_filepath ?? entry?.path);
      if (fromEntry) return fromEntry;
    }
  }
  const channel = channelName(model?.ion_class);
  if (!channel || !model?.modelDB_dir || !model?.mod_filename) return null;
  const filename = model.mod_filename.endsWith('.mod') ? model.mod_filename : `${model.mod_filename}.mod`;
  return [channel, `${model.modelDB_dir}_${filename}`];
}

function compactOriginalModel(model) {
  const {
    mod_filepath,
    mod_filename,
    unique_modelDB_mod_id,
    modelDB_dir,
    ICG,
    ion_class,
    ['Supermodel 1']: omnimodel1,
    ['Supermodel 2']: omnimodel2,
    Year,
    Authors,
  } = model;
  return {
    mod_filepath,
    mod_filename,
    unique_modelDB_mod_id,
    modelDB_dir,
    ICG,
    ion_class,
    'Supermodel 1': omnimodel1,
    'Supermodel 2': omnimodel2,
    Year,
    Authors,
  };
}

function compactIdenticalModel(model) {
  const { mod_filepath, mod_filename, unique_modelDB_mod_id, modelDB_dir, ICG, Year } = model;
  return { mod_filepath, mod_filename, unique_modelDB_mod_id, modelDB_dir, ICG, Year };
}

async function writeAtomically(file, value) {
  const temporary = `${file}.${process.pid}.tmp`;
  await assertWritableWithin(REPO_ROOT, file, 'Visualizer asset');
  await assertWritableWithin(REPO_ROOT, temporary, 'Temporary visualizer asset');
  await fs.mkdir(path.dirname(file), { recursive: true });
  let created = false;
  try {
    const handle = await fs.open(temporary, 'wx');
    created = true;
    try {
      await handle.writeFile(JSON.stringify(value));
    } finally {
      await handle.close();
    }
    await fs.rename(temporary, file);
  } finally {
    if (created) await fs.rm(temporary, { force: true });
  }
}

const network = JSON.parse(await fs.readFile(NETWORK_SOURCE, 'utf8'));
if (!Array.isArray(network?.nodes) || !Array.isArray(network?.links)) {
  throw new Error('network_data.json must contain node and link arrays');
}

const initialNodes = network.nodes
  .filter((node) => Number(node?.num_of_identicals ?? 0) >= INITIAL_MINIMUM_COPIES)
  .map(({ original_model, identical_models, num_of_identicals, label, id }) => ({
    original_model: compactOriginalModel(original_model),
    identical_models: identical_models.map(compactIdenticalModel),
    num_of_identicals,
    label,
    id,
    omnimodel_location: omnimodelLocation(original_model),
  }));
const initialNodeIds = new Set(initialNodes.map((node) => identifier(node.id)));
if (initialNodeIds.size !== initialNodes.length || initialNodes.length === 0) {
  throw new Error('Initial visualizer dataset has missing or duplicate node identifiers');
}

const initialLinks = network.links
  .filter((link) => initialNodeIds.has(identifier(link?.source)) && initialNodeIds.has(identifier(link?.target)))
  .map(({ source, target, weight }) => ({ source, target, weight }));
if (initialLinks.length === 0) {
  throw new Error('Initial visualizer dataset has no links');
}

await writeAtomically(INITIAL_NETWORK_OUTPUT, { nodes: initialNodes, links: initialLinks });

const layouts = JSON.parse(await fs.readFile(LAYOUT_SOURCE, 'utf8'));
for (let groupCount = MINIMUM_LAYOUT_GROUPS; groupCount <= MAXIMUM_LAYOUT_GROUPS; groupCount += 1) {
  const candidates = layouts?.[groupCount];
  if (
    !Array.isArray(candidates)
    || candidates.length === 0
    || candidates.some((candidate) => (
      !Array.isArray(candidate)
      || candidate.length !== groupCount
      || candidate.some((circle) => (
        !Array.isArray(circle)
        || circle.length !== 3
        || circle.some((value) => typeof value !== 'number' || !Number.isFinite(value))
      ))
    ))
  ) {
    throw new Error(`Invalid fixed-layout candidates for ${groupCount} groups`);
  }
  await writeAtomically(path.join(LAYOUT_OUTPUT_ROOT, `${groupCount}.json`), candidates);
}

const initialBytes = (await fs.stat(INITIAL_NETWORK_OUTPUT)).size;
console.log(
  `visualizer assets: ${initialNodes.length} initial nodes, ${initialLinks.length} initial links, `
  + `${Math.round(initialBytes / 1024)} KiB, ${MAXIMUM_LAYOUT_GROUPS - MINIMUM_LAYOUT_GROUPS + 1} layout tables`,
);

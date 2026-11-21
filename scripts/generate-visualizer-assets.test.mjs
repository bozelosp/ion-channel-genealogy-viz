import assert from 'node:assert/strict';
import { readFile } from 'node:fs/promises';
import test from 'node:test';
import { searchGraph } from '../lib/graphSearch.ts';

const readData = async (name) => JSON.parse(await readFile(new URL(`../public/data/${name}`, import.meta.url), 'utf8'));
const full = await readData('network_data.json');
const initial = await readData('network_data_initial.json');
const nodeIds = new Set(initial.nodes.map((node) => node.id));

test('filename searches include original and identical-copy filenames within the visible graph', () => {
  const nodes = [
    { id: 1, original_model: { mod_filename: 'H.mod' } },
    { id: 2, original_model: { mod_filename: 'ih.mod' }, identical_models: [{ mod_filename: 'na.mod' }] },
    { id: 3, original_model: { mod_filename: 'h.mod' } },
  ];
  const visible = new Set(['1', '2']);
  assert.deepEqual([...searchGraph(nodes, visible, ' h.MOD ')], ['1']);
  assert.deepEqual([...searchGraph(nodes, visible, 'na.mod')], ['2']);
  assert.equal(searchGraph(nodes, visible, 'na').size, 0);
  assert.equal(searchGraph(nodes, new Set(), 'h.mod').size, 0);
});

test('the full graph can find ModelDB entries outside the two-copy dataset', () => {
  const ids = (nodes) => new Set(nodes.map((node) => String(node.id)));
  const node = full.nodes.find((node) => node.original_model.unique_modelDB_mod_id === '119266_kaprox-ID-1');
  assert.ok(node);
  assert.ok(searchGraph(full.nodes, ids(full.nodes), '119266').has(String(node.id)));
  assert.equal(searchGraph(initial.nodes, ids(initial.nodes), '119266').has(String(node.id)), false);
});

test('initial graph retains every two-copy node and all edges between them, without changing weights', () => {
  const expected = full.nodes.filter((node) => (node.num_of_identicals || 1) >= 2);
  assert.deepEqual(initial.nodes.map((node) => node.id), expected.map((node) => node.id));
  assert.equal(nodeIds.size, expected.length);
  assert.deepEqual(initial.links, full.links
    .filter((link) => nodeIds.has(link.source) && nodeIds.has(link.target))
    .map(({ source, target, weight }) => ({ source, target, weight })));
});

test('compact metadata retains displayed, searchable, exported, and ancestor-selection fields', () => {
  const originals = new Map(full.nodes.map((node) => [node.id, node]));
  for (const node of initial.nodes) {
    const original = originals.get(node.id);
    assert.equal(node.label, original.label);
    assert.equal(node.num_of_identicals, original.num_of_identicals);
    assert.equal(node.identical_models.length, original.identical_models.length);
    for (const key of ['mod_filepath', 'mod_filename', 'unique_modelDB_mod_id', 'modelDB_dir', 'ICG', 'ion_class', 'Supermodel 1', 'Supermodel 2', 'Year', 'Authors']) {
      assert.deepEqual(node.original_model[key], original.original_model[key], `${node.id}: ${key}`);
    }
    node.identical_models.forEach((model, index) => {
      for (const key of ['mod_filepath', 'mod_filename', 'unique_modelDB_mod_id', 'modelDB_dir', 'ICG', 'Year']) {
        assert.deepEqual(model[key], original.identical_models[index][key], `${node.id}: copy ${index}: ${key}`);
      }
    });
  }
});

test('every ancestor candidate of an initial node is present in the initial graph', () => {
  const completeIds = new Set(full.nodes.map((node) => node.original_model.unique_modelDB_mod_id));
  const initialIds = new Set(initial.nodes.map((node) => node.original_model.unique_modelDB_mod_id));
  for (const node of initial.nodes) {
    for (const model of node.identical_models) {
      if (completeIds.has(model.unique_modelDB_mod_id)) {
        assert(initialIds.has(model.unique_modelDB_mod_id), `Missing ancestor for ${node.id}`);
      }
    }
  }
});

test('per-group layout tables retain every candidate and exact coordinate', async () => {
  const layouts = await readData('fixed_location_circles.json');
  for (let count = 2; count <= 25; count++) {
    assert.deepEqual(await readData(`fixed-location-circles/${count}.json`), layouts[count]);
  }
});

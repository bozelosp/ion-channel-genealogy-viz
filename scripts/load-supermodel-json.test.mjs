import assert from 'node:assert/strict';
import { spawnSync } from 'node:child_process';
import { promises as fs } from 'node:fs';
import os from 'node:os';
import path from 'node:path';
import { fileURLToPath } from 'node:url';
import test from 'node:test';
import {
  loadSupermodelJson,
  MAX_SUPERMODEL_JSON_BYTES,
  MAX_SUPERMODEL_JSON_CONTAINER_ENTRIES,
  MAX_SUPERMODEL_JSON_DEPTH,
  MAX_SUPERMODEL_JSON_VALUES,
} from './load-supermodel-json.mjs';

const GENERATOR = fileURLToPath(new URL('./generate-omnimodels.mjs', import.meta.url));
const RENDERER = fileURLToPath(new URL('./render-omnimodel-plots.mjs', import.meta.url));
// A protocol-4 pickle containing only an empty dictionary; no executable objects.
const EMPTY_DICTIONARY_PICKLE = Buffer.from([0x80, 0x04, 0x7d, 0x94, 0x2e]);

async function fixture(t) {
  const root = await fs.realpath(await fs.mkdtemp(path.join(os.tmpdir(), 'icg-supermodel-test-')));
  t.after(() => fs.rm(root, { recursive: true, force: true }));
  return {
    root,
    file: (...segments) => path.join(root, ...segments),
  };
}

function run(script, args, cwd) {
  return spawnSync(process.execPath, [script, ...args], {
    cwd,
    encoding: 'utf8',
  });
}

test('regular JSON remains usable through both generator entrypoints', async (t) => {
  const f = await fixture(t);
  const source = f.file('supermodel_data.json');
  const markdown = f.file('markdown');
  const reports = f.file('reports');
  const index = f.file('index.json');
  await fs.mkdir(markdown);
  await fs.writeFile(source, JSON.stringify({ Na: { Example: {} } }));

  assert.deepEqual(loadSupermodelJson(source), { Na: { Example: {} } });

  const generated = run(GENERATOR, [
    '--source', source,
    '--out-markdown', markdown,
    '--out-data', index,
    '--out-reports', reports,
  ], f.root);
  assert.equal(generated.status, 0, generated.stderr || generated.stdout);
  assert.deepEqual(JSON.parse(await fs.readFile(index, 'utf8')), [{
    channel: 'Na',
    originalSlug: 'Example',
    routeSlug: 'example',
    title: 'Example',
    summary: { suffix: null, gmaxName: null, states: [], gates: null },
  }]);

  const rendered = run(RENDERER, [
    '--source', source,
    '--out-markdown', markdown,
    '--dry-run',
  ], f.root);
  assert.equal(rendered.status, 0, rendered.stderr || rendered.stdout);
  assert.match(rendered.stdout, /No plots queued/);
});

test('binary pickle input is rejected before either caller writes output', async (t) => {
  const f = await fixture(t);
  const source = f.file('unsupported.pkl');
  const markdown = f.file('markdown');
  const reports = f.file('reports');
  const index = f.file('index.json');
  await fs.writeFile(source, EMPTY_DICTIONARY_PICKLE);

  assert.throws(() => loadSupermodelJson(source), /source is not valid (?:UTF-8|JSON)/);

  const generated = run(GENERATOR, [
    '--source', source,
    '--out-markdown', markdown,
    '--out-data', index,
    '--out-reports', reports,
  ], f.root);
  assert.equal(generated.status, 1, generated.stderr || generated.stdout);
  assert.match(generated.stderr, /source is not valid (?:UTF-8|JSON)/);
  await assert.rejects(fs.access(index));
  await assert.rejects(fs.access(reports));

  const rendered = run(RENDERER, [
    '--source', source,
    '--out-markdown', markdown,
    '--dry-run',
  ], f.root);
  assert.equal(rendered.status, 1, rendered.stderr || rendered.stdout);
  assert.match(rendered.stderr, /source is not valid (?:UTF-8|JSON)/);
  await assert.rejects(fs.access(markdown));
});

test('the shared loader refuses symlinks and files above the source-size bound', async (t) => {
  const f = await fixture(t);
  const source = f.file('source.json');
  const symlink = f.file('source-link.json');
  const oversized = f.file('oversized.json');
  await fs.writeFile(source, '{}');
  await fs.symlink(source, symlink);
  await fs.writeFile(oversized, '');
  await fs.truncate(oversized, MAX_SUPERMODEL_JSON_BYTES + 1);

  assert.throws(() => loadSupermodelJson(symlink), /regular, non-symlink/);
  assert.throws(() => loadSupermodelJson(oversized), /32 MiB limit/);
});

test('the shared loader rejects resource amplification before parsing', async (t) => {
  const f = await fixture(t);
  const tooDeep = f.file('too-deep.json');
  const tooManyEntries = f.file('too-many-entries.json');
  const tooManyValues = f.file('too-many-values.json');
  await fs.writeFile(
    tooDeep,
    '['.repeat(MAX_SUPERMODEL_JSON_DEPTH + 1) + '0' + ']'.repeat(MAX_SUPERMODEL_JSON_DEPTH + 1),
  );
  await fs.writeFile(
    tooManyEntries,
    `[${'[],'.repeat(MAX_SUPERMODEL_JSON_CONTAINER_ENTRIES)}[]]`,
  );
  const fullContainer = `[${'0,'.repeat(MAX_SUPERMODEL_JSON_CONTAINER_ENTRIES - 1)}0]`;
  const containerCount = Math.floor(MAX_SUPERMODEL_JSON_VALUES / MAX_SUPERMODEL_JSON_CONTAINER_ENTRIES) + 1;
  await fs.writeFile(
    tooManyValues,
    `[${Array.from({ length: containerCount }, () => fullContainer).join(',')}]`,
  );

  assert.throws(() => loadSupermodelJson(tooDeep), /nesting-depth limit/);
  assert.throws(() => loadSupermodelJson(tooManyEntries), /container limit/);
  assert.throws(() => loadSupermodelJson(tooManyValues), /value limit/);
});

test('the shared loader rejects malformed UTF-8 and non-finite JSON numbers', async (t) => {
  const f = await fixture(t);
  const malformedUtf8 = f.file('malformed-utf8.json');
  const nonFinite = f.file('non-finite.json');
  await fs.writeFile(malformedUtf8, Buffer.concat([
    Buffer.from('{"Na":{"Ex'),
    Buffer.from([0xc0]),
    Buffer.from('ample":{}}}'),
  ]));
  await fs.writeFile(nonFinite, '{"Na":{"Example":{"GMAX":1e400}}}');

  assert.throws(() => loadSupermodelJson(malformedUtf8), /not valid UTF-8/);
  assert.throws(() => loadSupermodelJson(nonFinite), /numbers must be finite/);
});

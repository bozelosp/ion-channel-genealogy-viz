import assert from 'node:assert/strict';
import { spawn, spawnSync } from 'node:child_process';
import { randomUUID } from 'node:crypto';
import { promises as fs } from 'node:fs';
import os from 'node:os';
import path from 'node:path';
import { fileURLToPath } from 'node:url';
import test from 'node:test';

// Exercise the real build coordinator, but never move this checkout's source
// or run Next.js. The tiny export/child fixtures make interruption tests cheap.
const SCRIPT = fileURLToPath(new URL('./build-static.mjs', import.meta.url));
const EXCLUDED = ['app/api/source-code', 'app/omnimodels/[channel]/[model]', 'app/omnimodels/static/[channel]/[slug]'];
const sleep = (ms) => new Promise((resolve) => setTimeout(resolve, ms));
const exists = (file) => fs.lstat(file).then(() => true, (error) => {
  if (error.code === 'ENOENT') return false;
  throw error;
});

async function until(check) {
  const deadline = Date.now() + 10_000;
  while (!(await check())) {
    if (Date.now() > deadline) throw new Error('Fixture did not reach expected state');
    await sleep(20);
  }
}

function alive(pid) {
  try { process.kill(pid, 0); return true; } catch (error) {
    if (error.code === 'ESRCH') return false;
    throw error;
  }
}

async function fixture(t) {
  const root = await fs.realpath(await fs.mkdtemp(path.join(os.tmpdir(), 'icg-static-test-')));
  const children = [];
  const file = (...segments) => path.join(root, ...segments);
  const write = async (relative, content) => {
    await fs.mkdir(path.dirname(file(relative)), { recursive: true });
    await fs.writeFile(file(relative), content);
  };
  t.after(async () => {
    for (const child of children) if (child.exitCode === null && child.signalCode === null) child.kill('SIGKILL');
    // Only this fixture's stub PID is ever targeted, including an orphan case.
    if (await exists(file('.stub-pid'))) {
      const pid = Number(await fs.readFile(file('.stub-pid'), 'utf8'));
      if (alive(pid)) process.kill(pid, 'SIGKILL');
    }
    await fs.rm(root, { recursive: true, force: true });
  });
  await write('package.json', '{"type":"module"}');
  await write('scripts/generate-shell-index.mjs', `
    import fs from 'node:fs';
    fs.writeFileSync('public/data/omnimodel-shell-index.json', JSON.stringify({count: 1}));
  `);
  await write('scripts/generate-visualizer-assets.mjs', `
    import fs from 'node:fs';
    fs.writeFileSync('public/data/network_data_initial.json', '{"nodes":[],"links":[]}');
    fs.mkdirSync('public/data/fixed-location-circles', {recursive: true});
    for (let count = 2; count <= 25; count++) {
      fs.writeFileSync('public/data/fixed-location-circles/' + count + '.json', '[]');
    }
  `);
  await fs.copyFile(SCRIPT, file('scripts/build-static.mjs'));
  await fs.copyFile(new URL('./run-build-child.mjs', import.meta.url), file('scripts/run-build-child.mjs'));
  await write('public/data/omnimodel-index.json', '[{"channel":"Na"}]');
  for (const relative of EXCLUDED) await write(`${relative}/route.ts`, 'original source\n');
  await write('node_modules/next/dist/bin/next', `
    import fs from 'node:fs';
    import path from 'node:path';
    const write = (name, data) => { fs.mkdirSync(path.dirname(name), {recursive: true}); fs.writeFileSync(name, data); };
    for (const name of ['index.html','404.html','robots.txt','data/omnimodel-index.json','data/network_data.json','data/network_data_initial.json','data/fixed_location_circles.json','data/network_data_16_oct_2019_1.json','data/omnimodel-reports/one.json','data/omnimodels/plot.png']) write('out/'+name, 'fixture');
    for (let count = 2; count <= 25; count++) write('out/data/fixed-location-circles/' + count + '.json', 'fixture');
    write('out/data/omnimodel-shell-index.json', '{"count":1}');
    if (process.env.STUB_HOLD === '1') {
      fs.writeFileSync('.stub-pid', String(process.pid));
      process.on('SIGTERM', () => {
        fs.writeFileSync('.stub-stopping', String(!fs.existsSync('app/api/source-code')));
        setTimeout(() => process.exit(0), 200);
      });
      setInterval(() => {}, 100);
    }
  `);
  function run(env = {}, nodeArgs = []) {
    const child = spawn(process.execPath, [...nodeArgs, file('scripts/build-static.mjs')], {
      cwd: root, env: { ...process.env, ...env }, stdio: ['ignore', 'pipe', 'pipe'],
    });
    children.push(child);
    let output = '';
    child.stdout.on('data', (data) => { output += data; });
    child.stderr.on('data', (data) => { output += data; });
    const result = new Promise((resolve, reject) => {
      child.on('error', reject);
      child.on('exit', (code, signal) => resolve({ code, signal, output }));
    });
    return { child, result };
  }
  async function intact(expected = 'original source\n') {
    for (const relative of EXCLUDED) assert.equal(await fs.readFile(file(relative, 'route.ts'), 'utf8'), expected);
    assert.equal(await exists(file('.static-excluded')), false);
    assert.equal(await exists(file('.static-build.lock')), false);
  }
  return { root, file, write, run, intact };
}

test('ordinary builds preserve uncommitted route edits and the canonical model count', async (t) => {
  const f = await fixture(t);
  assert.equal(spawnSync('git', ['init', '--quiet'], { cwd: f.root }).status, 0);
  assert.equal(spawnSync('git', ['add', 'app'], { cwd: f.root }).status, 0);
  for (const relative of EXCLUDED) await f.write(`${relative}/route.ts`, 'legitimate uncommitted edit\n');
  for (let run = 0; run < 2; run += 1) {
    const result = await f.run().result;
    assert.equal(result.code, 0, result.output);
    assert.match(result.output, /reports: 1, plots: 1, shell index count: 1/);
    assert.equal(await exists(f.file('out/data/fixed_location_circles.json')), true);
    assert.equal(await exists(f.file('out/data/network_data_16_oct_2019_1.json')), true);
    await f.intact('legitimate uncommitted edit\n');
  }
});

test('a rejected competitor cannot restore files or release the owner lock; SIGTERM waits for the child', async (t) => {
  const f = await fixture(t);
  const owner = f.run({ STUB_HOLD: '1' });
  await until(() => exists(f.file('.stub-pid')));
  const manifest = await fs.readFile(f.file('.static-excluded/manifest.json'), 'utf8');
  const lock = (await fs.readdir(f.file('.static-build.lock'))).sort();
  const competitor = await f.run().result;
  assert.equal(competitor.code, 1, competitor.output);
  assert.match(competitor.output, /another build/);
  assert.equal(await fs.readFile(f.file('.static-excluded/manifest.json'), 'utf8'), manifest);
  assert.deepEqual((await fs.readdir(f.file('.static-build.lock'))).sort(), lock);
  assert.equal(await exists(f.file(EXCLUDED[0])), false);
  owner.child.kill('SIGTERM');
  await until(() => exists(f.file('.stub-stopping')));
  assert.equal(await fs.readFile(f.file('.stub-stopping'), 'utf8'), 'true');
  assert.equal(await exists(f.file(EXCLUDED[0])), false);
  assert.equal(await exists(f.file('.static-build.lock')), true);
  const result = await owner.result;
  assert.equal(result.code, 143, result.output);
  await f.intact();
});

test('a crashed parent cannot be replaced while its child lives; a later run recovers', async (t) => {
  const f = await fixture(t);
  const owner = f.run({ STUB_HOLD: '1' });
  await until(() => exists(f.file('.stub-pid')));
  const childPid = Number(await fs.readFile(f.file('.stub-pid'), 'utf8'));
  owner.child.kill('SIGKILL');
  assert.equal((await owner.result).signal, 'SIGKILL');
  const competitor = await f.run().result;
  assert.equal(competitor.code, 1, competitor.output);
  assert.match(competitor.output, /another build/);
  assert.equal(await exists(f.file(EXCLUDED[0])), false);
  assert.equal(await fs.readFile(f.file('out/index.html'), 'utf8'), 'fixture');
  process.kill(childPid, 'SIGTERM');
  await until(() => !alive(childPid));
  const recovered = await f.run().result;
  assert.equal(recovered.code, 0, recovered.output);
  assert.match(recovered.output, /startup recovery/);
  await f.intact();
});

test('child failure restores sources', async (t) => {
  const f = await fixture(t);
  await f.write('scripts/generate-shell-index.mjs', 'process.exit(2);');
  const result = await f.run().result;
  assert.equal(result.code, 1, result.output);
  await f.intact();
});

test('a delayed stale reclaimer cannot delete a newly acquired live lock', async (t) => {
  const f = await fixture(t);
  const deadPid = spawnSync(process.execPath, ['-e', '']).pid;
  assert.equal(alive(deadPid), false);
  const staleName = `owner-${deadPid}-${randomUUID()}`;
  await f.write(`.static-build.lock/${staleName}`, '');
  // Deterministically pause at the real filesystem operation; no test hook
  // is shipped in the coordinator and no workspace source is moved.
  await f.write('pause-reclaimer.mjs', `
    import fs from 'node:fs';
    const unlink = fs.unlinkSync;
    fs.unlinkSync = function(file) {
      if (file === ${JSON.stringify(f.file('.static-build.lock', staleName))}) {
        fs.writeFileSync('.reclaimer-paused', '');
        const deadline = Date.now() + 10000;
        while (!fs.existsSync('.resume-reclaimer')) {
          if (Date.now() > deadline) throw new Error('Fixture resume timeout');
          Atomics.wait(new Int32Array(new SharedArrayBuffer(4)), 0, 0, 10);
        }
      }
      return unlink.call(fs, file);
    };
  `);
  const delayed = f.run({}, ['--import', f.file('pause-reclaimer.mjs')]);
  await Promise.race([
    until(() => exists(f.file('.reclaimer-paused'))),
    delayed.result.then((result) => { throw new Error(`Reclaimer exited before pause: ${result.output}`); }),
  ]);
  const owner = f.run({ STUB_HOLD: '1' });
  await until(() => exists(f.file('.stub-pid')));
  const liveRecords = (await fs.readdir(f.file('.static-build.lock'))).sort();
  await f.write('.resume-reclaimer', '');
  const loser = await delayed.result;
  assert.equal(loser.code, 1, loser.output);
  assert.deepEqual((await fs.readdir(f.file('.static-build.lock'))).sort(), liveRecords);
  assert.equal(await exists(f.file(EXCLUDED[0])), false);
  owner.child.kill('SIGTERM');
  assert.equal((await owner.result).code, 143);
  await f.intact();
});

test('a child whose unique owner vanished never loads build code', async (t) => {
  const f = await fixture(t);
  const result = spawnSync(process.execPath, [
    f.file('scripts/run-build-child.mjs'), f.file('.static-build.lock/missing-owner'),
    randomUUID(), f.file('node_modules/next/dist/bin/next'),
  ], { cwd: f.root, encoding: 'utf8' });
  assert.equal(result.status, 1);
  assert.match(result.stderr, /no longer owns the lock/);
  assert.equal(await exists(f.file('out')), false);
  await f.intact();
});

for (const mode of ['malformed', 'outside-root', 'unknown-file', 'source-conflict']) {
  test(`recovery preserves ${mode} state instead of deleting or overwriting data`, async (t) => {
    const f = await fixture(t);
    const target = f.file('.static-excluded/0-source-code');
    await fs.mkdir(path.dirname(target));
    const entry = { source: f.file(EXCLUDED[0]), target };
    if (mode === 'source-conflict') {
      await f.write('.static-excluded/0-source-code/route.ts', 'held source');
    } else if (mode === 'unknown-file') {
      await f.write('.static-excluded/unrelated.txt', 'preserve');
    } else if (mode === 'outside-root') {
      entry.source = path.join(f.root, '..', 'not-owned-by-this-build');
      await f.write('.static-excluded/0-source-code/route.ts', 'held source');
    }
    const manifest = mode === 'malformed' ? '[' : JSON.stringify([entry]);
    await f.write('.static-excluded/manifest.json', manifest);
    const result = await f.run().result;
    assert.equal(result.code, 1, result.output);
    assert.equal(await fs.readFile(f.file('.static-excluded/manifest.json'), 'utf8'), manifest);
    assert.equal(await fs.readFile(f.file(EXCLUDED[0], 'route.ts'), 'utf8'), 'original source\n');
    assert.equal(await exists(f.file('.static-build.lock')), false);
    if (mode === 'source-conflict' || mode === 'outside-root') assert.equal(await fs.readFile(path.join(target, 'route.ts'), 'utf8'), 'held source');
    if (mode === 'unknown-file') assert.equal(await fs.readFile(f.file('.static-excluded/unrelated.txt'), 'utf8'), 'preserve');
  });
}

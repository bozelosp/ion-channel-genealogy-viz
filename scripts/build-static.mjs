#!/usr/bin/env node
// Build the fully static production export.
//
// Static export cannot contain request-dependent server features, so the
// dev-only surfaces (the source-code proxy route handler and the two
// runtime-rendered report routes replaced by the client shell) are moved
// aside for the duration of the build and always restored afterwards.
//
// Interruption safety: the build takes an exclusive lock, records every move
// in .static-excluded/manifest.json BEFORE performing it, restores
// after its child exits on SIGINT/SIGTERM, recovers automatically
// from a previous run that died mid-build, and finishes by asserting that
// every excluded path is back. Existing working-tree edits are preserved.
// The export is then verified against the canonical tracked data.

import path from 'node:path';
import fsSync from 'node:fs';
import { promises as fs } from 'node:fs';
import { spawn } from 'node:child_process';
import { randomUUID } from 'node:crypto';
import { fileURLToPath } from 'node:url';

const REPO_ROOT = path.resolve(path.dirname(fileURLToPath(import.meta.url)), '..');
const HOLDING_DIR = path.join(REPO_ROOT, '.static-excluded');
const MANIFEST_PATH = path.join(HOLDING_DIR, 'manifest.json');
const LOCK_DIR = path.join(REPO_ROOT, '.static-build.lock');
const OUT_DIR = path.join(REPO_ROOT, 'out');

// Paths that only work on the dev server and must not exist during an export
// build: the dev source-code proxy route and the two runtime-rendered report
// routes that the exported client shell replaces in production.
const EXCLUDED_PATHS = [
  'app/api/source-code',
  'app/omnimodels/[channel]/[model]',
  'app/omnimodels/static/[channel]/[slug]',
];

function pathExistsSync(target) {
  try { fsSync.lstatSync(target); return true; } catch (error) {
    if (error.code === 'ENOENT') return false;
    throw error;
  }
}

// ---------------------------------------------------------------------------
// Synchronous moves and atomic manifest replacement keep signal handling
// from interleaving restoration with an in-flight source rename.

function readManifestSync() {
  if (!pathExistsSync(MANIFEST_PATH)) return [];
  const entries = JSON.parse(fsSync.readFileSync(MANIFEST_PATH, 'utf8'));
  const sources = new Set();
  if (!Array.isArray(entries) || entries.some((entry) => {
    const index = EXCLUDED_PATHS.findIndex((relative) => path.join(REPO_ROOT, relative) === entry?.source);
    if (index < 0 || entry.target !== heldPath(index) || sources.has(entry.source)) return true;
    sources.add(entry.source);
    return false;
  })) {
    throw new Error('Invalid exclusion manifest; preserved for manual recovery');
  }
  return entries;
}

function heldPath(index) {
  return path.join(HOLDING_DIR, `${index}-${path.basename(EXCLUDED_PATHS[index])}`);
}

function restoreSync(reason) {
  if (!ownsLock) return;
  if (!pathExistsSync(HOLDING_DIR)) return;
  const info = fsSync.lstatSync(HOLDING_DIR);
  if (!info.isDirectory() || info.isSymbolicLink()) {
    throw new Error('Exclusion state is not a real directory; preserved for manual recovery');
  }
  const entries = readManifestSync();
  const known = new Set(['manifest.json', 'manifest.json.next', ...entries.map((entry) => path.basename(entry.target))]);
  if (fsSync.readdirSync(HOLDING_DIR).some((name) => !known.has(name))) {
    throw new Error('Unknown exclusion files; preserved for manual recovery');
  }
  let failed = false;
  for (const { source, target } of [...entries].reverse()) {
    if (pathExistsSync(target)) {
      if (pathExistsSync(source)) {
        console.error(`restore conflict: both ${source} and ${target} exist - resolve manually`);
        failed = true;
        continue;
      }
      try {
        fsSync.renameSync(target, source);
      } catch (error) {
        console.error(`FAILED to restore ${source}: ${error.message}`);
        failed = true;
      }
    } else if (!pathExistsSync(source)) {
      console.error(`FAILED to restore ${source}: both paths are missing`);
      failed = true;
    }
  }
  if (!failed) {
    for (const name of ['manifest.json', 'manifest.json.next']) {
      const file = path.join(HOLDING_DIR, name);
      if (pathExistsSync(file)) fsSync.unlinkSync(file);
    }
    fsSync.rmdirSync(HOLDING_DIR);
    if (reason) console.error(`build-static: exclusions restored (${reason})`);
  } else {
    throw new Error('Exclusions could not be restored; preserved for manual recovery');
  }
}

let ownsLock = false;
const lockToken = randomUUID();
const ownerName = `owner-${process.pid}-${lockToken}`;
const ownerPath = path.join(LOCK_DIR, ownerName);

function releaseLockSync() {
  if (!ownsLock) return;
  ownsLock = false;
  removeRecordSync(ownerPath);
  try { fsSync.rmdirSync(LOCK_DIR); } catch (error) {
    // Another owner may have atomically installed a nonempty lock after we
    // removed our record. Never remove its files or its directory.
    if (!['ENOENT', 'ENOTEMPTY', 'EEXIST'].includes(error.code)) throw error;
  }
}

function removeRecordSync(file) {
  try { fsSync.unlinkSync(file); } catch (error) {
    if (error.code !== 'ENOENT') throw error;
  }
}

function processIsAlive(pid) {
  if (!Number.isSafeInteger(pid) || pid <= 1) throw new Error('Invalid build lock PID; leaving it untouched');
  try { process.kill(pid, 0); return true; } catch (error) {
    if (error.code === 'ESRCH') return false;
    throw error; // A permission failure is not proof that the owner is dead.
  }
}

// ---------------------------------------------------------------------------

function acquireLock() {
  // Publish a complete, nonempty directory in one atomic rename. Unlike
  // mkdir + a shared PID filename, no empty acquisition window or reused
  // record name lets a stale reclaimer delete a new owner's lock.
  const prepared = fsSync.mkdtempSync(`${LOCK_DIR}-`);
  const record = path.join(prepared, ownerName);
  fsSync.writeFileSync(record, '', { flag: 'wx', mode: 0o600 });
  try {
    try {
      fsSync.renameSync(prepared, LOCK_DIR);
    } catch (error) {
      if (!['EEXIST', 'ENOTEMPTY'].includes(error.code)) throw error;
      const info = fsSync.lstatSync(LOCK_DIR);
      if (!info.isDirectory() || info.isSymbolicLink()) throw new Error('Build lock is not a real directory');
      const names = fsSync.readdirSync(LOCK_DIR);
      for (const name of names) {
        const match = /^(?:owner|child)-(\d+)-[0-9a-f-]{36}$/.exec(name);
        const fileInfo = fsSync.lstatSync(path.join(LOCK_DIR, name));
        if (!match || !fileInfo.isFile() || fileInfo.isSymbolicLink()) {
          throw new Error('Invalid build lock record; leaving it untouched');
        }
        if (processIsAlive(Number(match[1]))) throw new Error(`another build (pid ${match[1]}) holds ${LOCK_DIR}`);
      }
      console.error('build-static: clearing stale lock from a dead build');
      // Only delete the exact dead records observed above. Every owner and
      // child has a unique filename; rmdir cannot remove a new nonempty lock.
      for (const name of names) removeRecordSync(path.join(LOCK_DIR, name));
      try { fsSync.rmdirSync(LOCK_DIR); } catch (error) {
        if (error.code !== 'ENOENT') throw error;
      }
      fsSync.renameSync(prepared, LOCK_DIR);
    }
    ownsLock = true;
  } finally {
    if (pathExistsSync(prepared)) {
      removeRecordSync(record);
      fsSync.rmdirSync(prepared);
    }
  }
}

// Children run via async spawn so the event loop stays free and signal
// handlers fire immediately (spawnSync would defer them until the child
// exited). Wait for the child to close before restoring its input paths.
let currentChild = null;
let interrupted = null;

function throwIfInterrupted() {
  if (interrupted) throw new Error(`Build interrupted by ${interrupted}`);
}

function runChild(args, env) {
  throwIfInterrupted();
  return new Promise((resolve, reject) => {
    currentChild = spawn(process.execPath, [
      path.join(REPO_ROOT, 'scripts', 'run-build-child.mjs'), ownerPath, lockToken, ...args,
    ], {
      cwd: REPO_ROOT,
      stdio: 'inherit',
      env: env ? { ...process.env, ...env } : process.env,
    });
    const childRecord = currentChild.pid ? path.join(LOCK_DIR, `child-${currentChild.pid}-${lockToken}`) : null;
    currentChild.on('error', reject);
    currentChild.on('close', (code, signal) => {
      currentChild = null;
      try {
        if (childRecord) removeRecordSync(childRecord);
        resolve(signal ? 1 : code);
      } catch (error) {
        reject(error);
      }
    });
  });
}

function installSafetyHandlers() {
  const bail = (signal) => {
    interrupted ??= signal;
    console.error(`\nbuild-static: received ${signal}, stopping build before restore`);
    if (currentChild) {
      try { currentChild.kill('SIGTERM'); } catch {}
    }
  };
  process.on('SIGINT', () => bail('SIGINT'));
  process.on('SIGTERM', () => bail('SIGTERM'));
}

async function run() {
  acquireLock();
  installSafetyHandlers();

  // Recover from a previous run that died with exclusions still held.
  if (pathExistsSync(HOLDING_DIR)) {
    console.error('build-static: found abandoned exclusion state, recovering');
    restoreSync('startup recovery');
    if (pathExistsSync(HOLDING_DIR)) {
      throw new Error(`${HOLDING_DIR} could not be auto-recovered - resolve manually`);
    }
  }

  fsSync.mkdirSync(HOLDING_DIR);
  const moved = [];

  for (const [index, relative] of EXCLUDED_PATHS.entries()) {
    const source = path.join(REPO_ROOT, relative);
    if (!pathExistsSync(source)) continue;
    const target = heldPath(index);
    // Record the intended move BEFORE performing it so recovery always
    // knows the mapping even if we die between rename and the next write.
    moved.push({ source, target });
    fsSync.writeFileSync(`${MANIFEST_PATH}.next`, JSON.stringify(moved));
    fsSync.renameSync(`${MANIFEST_PATH}.next`, MANIFEST_PATH);
    fsSync.renameSync(source, target);
  }
  console.log(`build-static: excluded ${moved.length} dev-only paths`);

  const shellIndex = await runChild([path.join(REPO_ROOT, 'scripts', 'generate-shell-index.mjs')]);
  throwIfInterrupted();
  if (shellIndex !== 0) throw new Error('shell index generation failed');
  const visualizerAssets = await runChild([path.join(REPO_ROOT, 'scripts', 'generate-visualizer-assets.mjs')]);
  throwIfInterrupted();
  if (visualizerAssets !== 0) throw new Error('visualizer asset generation failed');
  const expectedReportCount = JSON.parse(
    await fs.readFile(path.join(REPO_ROOT, 'public', 'data', 'omnimodel-index.json'), 'utf8'),
  ).length;

  // Stale generated route types from dev/server builds reference the
  // excluded paths, so export builds always start from a clean build dir.
  await fs.rm(path.join(REPO_ROOT, '.next'), { recursive: true, force: true });
  await fs.rm(OUT_DIR, { recursive: true, force: true });

  const build = await runChild(
    [path.join(REPO_ROOT, 'node_modules', 'next', 'dist', 'bin', 'next'), 'build', '--webpack'],
    { STATIC_EXPORT: '1' },
  );
  throwIfInterrupted();
  if (build !== 0) throw new Error('next build failed');

  restoreSync('build complete');

  // --- Post-restore assertions --------------------------------------------
  const unrestored = EXCLUDED_PATHS.filter(
    (relative) => moved.some((m) => m.source === path.join(REPO_ROOT, relative))
      && !pathExistsSync(path.join(REPO_ROOT, relative)),
  );
  if (unrestored.length > 0 || pathExistsSync(HOLDING_DIR)) {
    throw new Error(`exclusions not fully restored: ${unrestored.join(', ') || HOLDING_DIR}`);
  }
  // --- Post-process and verify the export ---------------------------------

  // Keep canonical scientific downloads available. The current application
  // requests only its compact graph and small group-count layout tables.

  const requiredFiles = [
    'index.html',
    '404.html',
    'robots.txt',
    'data/omnimodel-index.json',
    'data/omnimodel-shell-index.json',
    'data/network_data.json',
    'data/network_data_initial.json',
    'data/fixed_location_circles.json',
  ];
  for (const relative of requiredFiles) {
    if (!pathExistsSync(path.join(OUT_DIR, relative))) {
      throw new Error(`Export is missing required file: ${relative}`);
    }
  }
  for (let groupCount = 2; groupCount <= 25; groupCount += 1) {
    const relative = `data/fixed-location-circles/${groupCount}.json`;
    if (!pathExistsSync(path.join(OUT_DIR, relative))) {
      throw new Error(`Export is missing required file: ${relative}`);
    }
  }

  async function countFiles(dir, predicate) {
    let count = 0;
    const entries = await fs.readdir(dir, { withFileTypes: true });
    for (const entry of entries) {
      const full = path.join(dir, entry.name);
      if (entry.isDirectory()) count += await countFiles(full, predicate);
      else if (predicate(entry.name)) count += 1;
    }
    return count;
  }

  const reportCount = await countFiles(
    path.join(OUT_DIR, 'data', 'omnimodel-reports'),
    (name) => name.endsWith('.json'),
  );
  if (reportCount !== expectedReportCount) {
    throw new Error(`Expected ${expectedReportCount} report objects in export, found ${reportCount}`);
  }

  const shellIndexPayload = JSON.parse(
    await fs.readFile(path.join(OUT_DIR, 'data', 'omnimodel-shell-index.json'), 'utf8'),
  );
  if (shellIndexPayload.count !== expectedReportCount) {
    throw new Error(`Shell index count ${shellIndexPayload.count} != ${expectedReportCount}`);
  }

  const plotCount = await countFiles(
    path.join(OUT_DIR, 'data', 'omnimodels'),
    (name) => /\.(?:png|jpe?g|gif|svg|webp)$/i.test(name),
  );

  const htmlFiles = [];
  async function listHtml(dir) {
    const entries = await fs.readdir(dir, { withFileTypes: true });
    for (const entry of entries) {
      const full = path.join(dir, entry.name);
      if (entry.isDirectory()) {
        if (entry.name === 'data') continue;
        await listHtml(full);
      } else if (entry.name.endsWith('.html')) {
        htmlFiles.push(path.relative(OUT_DIR, full));
      }
    }
  }
  await listHtml(OUT_DIR);

  console.log('build-static: export verified');
  console.log(`  reports: ${reportCount}, plots: ${plotCount}, shell index count: ${shellIndexPayload.count}`);
  console.log(`  html pages: ${htmlFiles.sort().join(', ')}`);
}

try {
  await run();
} catch (error) {
  console.error(`build-static: ${error.message}`);
  process.exitCode = 1;
} finally {
  try {
    restoreSync(null);
  } finally {
    releaseLockSync();
  }
  if (interrupted) process.exitCode = interrupted === 'SIGINT' ? 130 : 143;
}

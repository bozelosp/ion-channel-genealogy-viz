#!/usr/bin/env node

// Synchronise generated Omnimodel markdown + plots into the icgenealogy
// channel repositories using the naming convention:
//   omnimodel.md
//   omni-plots/

import { promises as fs } from 'node:fs';
import path from 'node:path';
import { fileURLToPath } from 'node:url';
import { randomUUID } from 'node:crypto';
import { resolveWithin, safePathSegment } from './path-safety.mjs';

const __filename = fileURLToPath(import.meta.url);
const __dirname = path.dirname(__filename);
const projectRoot = path.resolve(__dirname, '..');
const appRoot = projectRoot;

function usage() {
  return `Usage: node scripts/sync-omnimodel-assets.mjs [options]\n\n` +
    `Options:\n` +
    `  --source <path>   Source root (default: public/data/omnimodels)\n` +
    `  --target <path>   Destination repo root (required)\n` +
    `  --channel <name>  Only process the given ion class (e.g. Ca, Na)\n` +
    `  --summary         Suppress per-model logs; show channel summary only\n` +
    `  --dry-run         Preview operations without copying\n` +
    `  --help            Show this message\n`;
}

function parseArgs() {
  const args = process.argv.slice(2);
  const opts = {
    sourceRoot: path.join(appRoot, 'public', 'data', 'omnimodels'),
    targetRoot: null,
    channelFilter: null,
    summaryOnly: false,
    dryRun: false,
  };

  for (let i = 0; i < args.length; i += 1) {
    const arg = args[i];
    if (arg === '--source' && args[i + 1]) {
      opts.sourceRoot = path.resolve(args[++i]);
    } else if (arg === '--target' && args[i + 1]) {
      opts.targetRoot = path.resolve(args[++i]);
    } else if (arg === '--channel' && args[i + 1]) {
      opts.channelFilter = args[++i];
    } else if (arg === '--summary') {
      opts.summaryOnly = true;
    } else if (arg === '--dry-run') {
      opts.dryRun = true;
    } else if (arg === '--help' || arg === '-h') {
      console.log(usage());
      process.exit(0);
    } else {
      throw new Error(`Unknown argument: ${arg}`);
    }
  }

  if (!opts.targetRoot) {
    console.error('Missing --target <path> pointing at the channel repo checkout.');
    process.exit(1);
  }
  return opts;
}

async function exists(p) {
  try {
    await fs.access(p);
    return true;
  } catch {
    return false;
  }
}

async function requireDirectory(directoryPath, label) {
  const info = await fs.lstat(directoryPath);
  if (!info.isDirectory() || info.isSymbolicLink()) {
    throw new Error(`${label} must be a real directory, not a symlink`);
  }
}

async function requireRegularFile(filePath, label) {
  const info = await fs.lstat(filePath);
  if (!info.isFile() || info.isSymbolicLink()) {
    throw new Error(`${label} must be a regular file, not a symlink`);
  }
}

async function copyFileAtomically(src, dest) {
  const temporary = path.join(
    path.dirname(dest),
    `.${path.basename(dest)}.${process.pid}.${randomUUID()}.tmp`,
  );
  try {
    await fs.copyFile(src, temporary);
    await fs.rename(temporary, dest);
  } finally {
    await fs.unlink(temporary).catch(() => {});
  }
}

async function copyMarkdown(src, dest, dryRun) {
  await requireRegularFile(src, 'Markdown source');
  if (dryRun) return;
  await copyFileAtomically(src, dest);
}

async function writePlots(srcDir, destDir, dryRun) {
  if (!(await exists(srcDir))) {
    return { files: 0 };
  }
  await requireDirectory(srcDir, 'Plot source');
  const sourceEntries = await fs.readdir(srcDir, { withFileTypes: true });
  const invalidSourceEntries = sourceEntries.filter(
    (entry) => !entry.isFile() || entry.isSymbolicLink() || !/\.png$/i.test(entry.name),
  );
  if (invalidSourceEntries.length > 0) {
    throw new Error(`Plot source contains unsupported entries: ${invalidSourceEntries.map((entry) => entry.name).join(', ')}`);
  }
  const files = sourceEntries.map((entry) => safePathSegment(entry.name, 'plot filename'));
  if (dryRun) {
    return { files: files.length };
  }

  if (await exists(destDir)) {
    await requireDirectory(destDir, 'Plot destination');
  } else {
    await fs.mkdir(destDir, { recursive: true });
  }

  const existingEntries = await fs.readdir(destDir, { withFileTypes: true });
  const managedExisting = existingEntries.filter((entry) => /\.png$/i.test(entry.name));
  const unsafeExisting = managedExisting.filter((entry) => !entry.isFile() || entry.isSymbolicLink());
  if (unsafeExisting.length > 0) {
    throw new Error(`Plot destination contains unsafe managed entries: ${unsafeExisting.map((entry) => entry.name).join(', ')}`);
  }

  await Promise.all(files.map((file) => copyFileAtomically(
    resolveWithin(srcDir, file),
    resolveWithin(destDir, file),
  )));

  const sourceNames = new Set(files);
  const staleFiles = managedExisting
    .map((entry) => safePathSegment(entry.name, 'existing plot filename'))
    .filter((file) => !sourceNames.has(file));
  await Promise.all(staleFiles.map((file) => fs.unlink(resolveWithin(destDir, file))));
  return { files: files.length };
}

async function main() {
  const opts = parseArgs();
  const { sourceRoot, targetRoot, channelFilter, summaryOnly, dryRun } = opts;

  await requireDirectory(sourceRoot, 'Source root');
  await requireDirectory(targetRoot, 'Target root');
  const safeChannelFilter = channelFilter ? safePathSegment(channelFilter, 'channel filter') : null;

  const channelEntries = await fs.readdir(sourceRoot, { withFileTypes: true });
  const summary = [];

  for (const channelDir of channelEntries) {
    if (!channelDir.isDirectory()) continue;
    const channel = safePathSegment(channelDir.name, 'channel');
    if (safeChannelFilter && channel !== safeChannelFilter) continue;

    const sourceChannel = resolveWithin(sourceRoot, channel);
    const channelPath = resolveWithin(targetRoot, channel);
    let targetChannel;
    if (await exists(channelPath)) {
      await requireDirectory(channelPath, 'Target channel');
      targetChannel = channelPath;
    } else if (safeChannelFilter === channel) {
      targetChannel = targetRoot;
    } else {
      throw new Error(`Target has no ${channel} directory; use --channel ${channel} for a single-channel repository`);
    }

    await requireDirectory(targetChannel, 'Target channel root');

    const modelEntries = await fs.readdir(sourceChannel, { withFileTypes: true });
    let synced = 0;
    let skipped = 0;
    const missing = [];

    for (const modelDir of modelEntries) {
      if (!modelDir.isDirectory()) continue;
      const slug = safePathSegment(modelDir.name, 'model slug');
      const sourceSlug = resolveWithin(sourceChannel, slug);
      const targetSlug = resolveWithin(targetChannel, slug);

      if (!(await exists(targetSlug))) {
        missing.push(slug);
        skipped += 1;
        continue;
      }
      await requireDirectory(sourceSlug, 'Source model directory');
      await requireDirectory(targetSlug, 'Target model directory');

      if (!summaryOnly) {
        console.log(`Syncing ${channel}/${slug}`);
      }

      const markdownSrc = resolveWithin(sourceSlug, `${slug}.md`);
      const markdownDest = resolveWithin(targetSlug, 'omnimodel.md');
      await copyMarkdown(markdownSrc, markdownDest, dryRun);

      const plotsSrc = resolveWithin(sourceSlug, 'images');
      const plotsDest = resolveWithin(targetSlug, 'omni-plots');
      const { files } = await writePlots(plotsSrc, plotsDest, dryRun);

      if (!summaryOnly && dryRun) {
        console.log(`  markdown -> ${markdownDest}`);
        if (files) {
          console.log(`  plots -> ${plotsDest} (${files} file${files === 1 ? '' : 's'})`);
        }
      }

      synced += 1;
    }

    summary.push({ channel, synced, skipped, missing });
  }

  console.log('Sync completed successfully.\n');
  for (const { channel, synced, skipped, missing } of summary) {
    const parts = [`Channel ${channel}: ${synced} synced`];
    if (skipped) parts.push(`${skipped} skipped`);
    console.log(parts.join(', '));
    if (missing.length > 0) {
      console.log(`  Missing folders (not updated): ${missing.join(', ')}`);
    }
  }
}

main().catch((error) => {
  console.error('Sync failed:', error.message ?? error);
  process.exit(1);
});

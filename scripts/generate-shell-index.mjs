#!/usr/bin/env node
// Generate public/data/omnimodel-shell-index.json: the compact lookup the static
// report shell fetches to resolve a model URL into its report JSON and plot
// assets without prerendering one page per model.
//
// Keys are lowercased "<channel>/<originalSlug>"; byRoute maps the
// "/omnimodels/static/<channel>/<routeSlug>" address space onto the same
// entries. Values stay positional arrays to keep the payload small.

import path from 'node:path';
import { promises as fs } from 'node:fs';
import { fileURLToPath } from 'node:url';

const REPO_ROOT = path.resolve(path.dirname(fileURLToPath(import.meta.url)), '..');
const INDEX_PATH = path.join(REPO_ROOT, 'public', 'data', 'omnimodel-index.json');
const ASSETS_ROOT = path.join(REPO_ROOT, 'public', 'data', 'omnimodels');
const OUTPUT_PATH = path.join(REPO_ROOT, 'public', 'data', 'omnimodel-shell-index.json');
const IMAGE_EXTENSION = /\.(?:png|jpe?g|gif|svg|webp)$/i;

function shellKey(channel, slug) {
  return `${channel.toLowerCase()}/${slug.toLowerCase()}`;
}

async function listImages(channel, originalSlug) {
  const imagesDir = path.join(ASSETS_ROOT, channel, originalSlug, 'images');
  try {
    const entries = await fs.readdir(imagesDir, { withFileTypes: true });
    return entries
      .filter((entry) => entry.isFile() && IMAGE_EXTENSION.test(entry.name))
      .map((entry) => entry.name)
      .sort((a, b) => a.localeCompare(b));
  } catch (error) {
    if (error && error.code === 'ENOENT') return [];
    throw error;
  }
}

const index = JSON.parse(await fs.readFile(INDEX_PATH, 'utf8'));
if (!Array.isArray(index) || index.length === 0) {
  throw new Error('omnimodel-index.json is empty or not an array');
}

const models = {};
const byRoute = {};

for (const item of index) {
  const { channel, originalSlug, routeSlug, title } = item;
  if (!channel || !originalSlug || !routeSlug || !title) {
    throw new Error(`Incomplete omnimodel index entry: ${JSON.stringify(item)}`);
  }
  const key = shellKey(channel, originalSlug);
  if (models[key]) throw new Error(`Duplicate shell key: ${key}`);
  models[key] = [channel, originalSlug, routeSlug, title, await listImages(channel, originalSlug)];

  const routeKey = shellKey(channel, routeSlug);
  if (byRoute[routeKey] && byRoute[routeKey] !== key) {
    throw new Error(`Route key collision: ${routeKey}`);
  }
  byRoute[routeKey] = key;
}

const payload = { count: index.length, models, byRoute };
await fs.writeFile(OUTPUT_PATH, JSON.stringify(payload));

const written = await fs.stat(OUTPUT_PATH);
console.log(`shell index: ${index.length} models, ${Object.keys(byRoute).length} route aliases, ${(written.size / 1024).toFixed(0)} KiB`);

#!/usr/bin/env node

import { spawnSync } from 'node:child_process';
import { promises as fs } from 'node:fs';
import { randomUUID } from 'node:crypto';
import path from 'node:path';
import { fileURLToPath } from 'node:url';
import {
  assertRemovableWithin,
  assertWritableWithin,
  ensureRealDirectory,
  resolveWithin,
  safePathSegment,
} from './path-safety.mjs';

const __filename = fileURLToPath(import.meta.url);
const __dirname = path.dirname(__filename);
const projectRoot = path.resolve(__dirname, '..');

const DEFAULT_SOURCE = path.resolve(projectRoot, 'supermodels', 'supermodel_data.pkl');
const DEFAULT_MARKDOWN_ROOT = path.join(projectRoot, 'public', 'data', 'omnimodels');
const DEFAULT_GENERATED_DATA = path.join(projectRoot, 'public', 'data', 'omnimodel-index.json');
const DEFAULT_REPORTS_ROOT = path.join(projectRoot, 'public', 'data', 'omnimodel-reports');
const GENERATED_DATA_ROOT = path.join(projectRoot, 'public', 'data');
const PICKLE_LOADER = path.join(__dirname, 'load-supermodel-json.py');

function formatDuration(seconds) {
  if (!Number.isFinite(seconds)) return '?:??';
  const s = Math.max(0, Math.round(seconds));
  const m = Math.floor(s / 60);
  const rem = s % 60;
  if (m >= 60) {
    const h = Math.floor(m / 60);
    const mm = m % 60;
    return `${h}h ${mm}m`;
  }
  return `${m}:${rem.toString().padStart(2, '0')}`;
}

function createProgressTracker(total, label) {
  const start = Date.now();
  function update(current) {
    if (total === 0) return;
    const elapsed = (Date.now() - start) / 1000;
    const rate = current / Math.max(elapsed, 0.001);
    const remaining = total - current;
    const eta = remaining / Math.max(rate, 0.001);
    const percent = ((current / total) * 100).toFixed(1);
    const message = `${label}: ${current}/${total} (${percent}%) | ETA ${formatDuration(eta)}`;
    process.stdout.write(`\r${message.padEnd(80, ' ')}`);
    if (current === total) {
      process.stdout.write('\n');
    }
  }
  return { update };
}

function parseArgs() {
  const options = {
    source: DEFAULT_SOURCE,
    markdownDir: DEFAULT_MARKDOWN_ROOT,
    generatedDataFile: DEFAULT_GENERATED_DATA,
    reportsRoot: DEFAULT_REPORTS_ROOT,
    cleanData: false,
  };

  const args = process.argv.slice(2);
  for (let i = 0; i < args.length; i += 1) {
    const arg = args[i];
    if (arg === '--source' && args[i + 1]) {
      options.source = path.resolve(process.cwd(), args[i + 1]);
      i += 1;
    } else if (arg === '--out-markdown' && args[i + 1]) {
      options.markdownDir = path.resolve(process.cwd(), args[i + 1]);
      i += 1;
    } else if (arg === '--out-data' && args[i + 1]) {
      options.generatedDataFile = path.resolve(process.cwd(), args[i + 1]);
      i += 1;
    } else if (arg === '--out-reports' && args[i + 1]) {
      options.reportsRoot = path.resolve(process.cwd(), args[i + 1]);
      i += 1;
    } else if (arg === '--clean') {
      options.cleanData = true;
    } else if (arg === '--help' || arg === '-h') {
      printUsage();
      process.exit(0);
    } else {
      throw new Error(`Unknown argument: ${arg}`);
    }
  }

  return options;
}

function printUsage() {
  console.log(`Usage: node scripts/generate-omnimodels.mjs [options]

Options:
  --source <path>         Path to supermodel_data.pkl (default: ${DEFAULT_SOURCE})
  --out-markdown <path>   Directory to write markdown output (default: ${DEFAULT_MARKDOWN_ROOT})
  --out-data <path>       Path to compact generated index JSON (default: ${DEFAULT_GENERATED_DATA})
  --out-reports <path>    Directory for per-model report JSON (default: ${DEFAULT_REPORTS_ROOT})
  --clean                 Remove previously generated data file before writing a new one
  -h, --help              Show this help message
`);
}

function runPythonExtraction(sourcePath) {
  const result = spawnSync('python3', [PICKLE_LOADER, sourcePath], {
    encoding: 'utf-8',
    maxBuffer: 128 * 1024 * 1024,
    timeout: 60_000,
  });

  if (result.error) {
    throw result.error;
  }

  if (result.status !== 0) {
    throw new Error(`Python exited with status ${result.status}: ${result.stderr}`);
  }

  return JSON.parse(result.stdout);
}

function formatHeading(key) {
  return key.replace(/_/g, ' ');
}

function formatValueAsMarkdown(value, indentLevel = 0) {
  const indent = ' '.repeat(indentLevel);
  if (value === null || value === undefined) {
    return `${indent}null`;
  }

  const valueType = typeof value;
  if (valueType === 'string') {
    return `${indent}${value}`;
  }

  if (valueType === 'number' || valueType === 'boolean') {
    return `${indent}${String(value)}`;
  }

  if (Array.isArray(value)) {
    if (value.length === 0) {
      return `${indent}[]`;
    }
    const simple = value.every((item) => {
      const type = typeof item;
      return item === null || type === 'number' || type === 'string' || type === 'boolean';
    });
    if (simple && value.length <= 20) {
      return value.map((item) => `${indent}- ${item}`).join('\n');
    }
    const jsonBlock = JSON.stringify(value, null, 2);
    return `${indent}\`\`\`json\n${jsonBlock}\n${indent}\`\`\``;
  }

  if (valueType === 'object') {
    const entries = Object.entries(value);
    if (entries.length === 0) {
      return `${indent}{}`;
    }
    const simple = entries.every(([, item]) => {
      const type = typeof item;
      return item === null || type === 'number' || type === 'string' || type === 'boolean';
    });
    if (simple && entries.length <= 20) {
      return entries
        .map(([k, v]) => `${indent}- **${k}**: ${v === null ? 'null' : v}`)
        .join('\n');
    }
    const jsonBlock = JSON.stringify(value, null, 2);
    return `${indent}\`\`\`json\n${jsonBlock}\n${indent}\`\`\``;
  }

  return `${indent}${String(value)}`;
}

async function ensureDir(dirPath) {
  await fs.mkdir(dirPath, { recursive: true });
}

async function writeFileAtomically(root, targetPath, contents) {
  const safeTarget = await assertWritableWithin(root, targetPath, 'generated output');
  const temporaryPath = path.join(
    path.dirname(safeTarget),
    `.${path.basename(safeTarget)}.${process.pid}.${randomUUID()}.tmp`,
  );
  try {
    await fs.writeFile(temporaryPath, contents, { flag: 'wx', mode: 0o644 });
    await fs.rename(temporaryPath, safeTarget);
  } finally {
    await fs.unlink(temporaryPath).catch(() => {});
  }
}

function sanitizeRouteSlug(raw) {
  const cleaned = raw
    .toLowerCase()
    .replace(/[^a-z0-9]+/g, '-')
    .replace(/^-+|-+$/g, '');
  return cleaned || 'omnimodel';
}

function ensureUniqueSlug(map, channel, slug) {
  const key = channel.toLowerCase();
  if (!map.has(key)) {
    map.set(key, new Set());
  }
  const slugSet = map.get(key);
  let candidate = slug;
  let counter = 1;
  while (slugSet.has(candidate)) {
    counter += 1;
    candidate = `${slug}-${counter}`;
  }
  slugSet.add(candidate);
  return candidate;
}

function toStringArray(value) {
  if (!Array.isArray(value)) return [];
  return value.map((item) => String(item));
}

function toNumericRecord(value) {
  if (!value || typeof value !== 'object' || Array.isArray(value)) {
    return null;
  }
  const result = {};
  for (const [key, val] of Object.entries(value)) {
    const numeric = typeof val === 'number' ? val : Number(val);
    result[key] = Number.isFinite(numeric) ? numeric : null;
  }
  return result;
}

async function collectImageMarkdown(markdownDir) {
  const imagesDir = path.join(markdownDir, 'images');
  try {
    const entries = await fs.readdir(imagesDir, { withFileTypes: true });
    const imageFiles = entries
      .filter((entry) => entry.isFile())
      .map((entry) => entry.name)
      .filter((name) => /\.(png|jpg|jpeg|gif|svg)$/i.test(name))
      .sort((a, b) => a.localeCompare(b));
    if (imageFiles.length === 0) {
      return '';
    }
    const lines = ['## Figures', ''];
    for (const file of imageFiles) {
      lines.push(`![${file}](images/${file})`);
    }
    lines.push('');
    return lines.join('\n');
  } catch (error) {
    if (error && error.code !== 'ENOENT') {
      console.warn(`Failed to read images for ${markdownDir}:`, error.message);
    }
    return '';
  }
}

function buildMarkdown(channel, modelName, data, imageSection) {
  const lines = [`# ${modelName}`, '', `**Ion Channel Class**: ${channel}`, ''];

  for (const [key, value] of Object.entries(data)) {
    lines.push(`## ${formatHeading(key)}`);
    lines.push('');
    lines.push(formatValueAsMarkdown(value));
    lines.push('');
  }

  if (imageSection) {
    lines.push(imageSection.trimEnd());
    lines.push('');
  }

  return lines.join('\n');
}

function buildReport(channel, modelName, routeSlug, data) {
  const sections = Object.entries(data).map(([key, value]) => ({
    key,
    title: formatHeading(key),
    value,
  }));

  return {
    channel,
    originalSlug: modelName,
    routeSlug,
    title: modelName,
    summary: {
      suffix: data.SUFFIX ?? null,
      gmaxName: data.GMAX_NAME ?? null,
      states: toStringArray(data.STATES),
      gates: toNumericRecord(data.GATES),
    },
    sections,
  };
}

async function writeMarkdownFile(markdownRoot, channel, modelName, markdownContent) {
  const safeChannel = safePathSegment(channel, 'channel');
  const safeModel = safePathSegment(modelName, 'model name');
  const targetDir = resolveWithin(markdownRoot, safeChannel, safeModel);
  await assertWritableWithin(markdownRoot, targetDir, 'markdown output directory');
  await ensureDir(targetDir);
  const markdownPath = resolveWithin(targetDir, `${safeModel}.md`);
  await writeFileAtomically(markdownRoot, markdownPath, markdownContent);
  return markdownPath;
}

async function writeGeneratedDataFile(generatedDataFile, reports) {
  const outputRoot = await ensureRealDirectory(path.dirname(generatedDataFile), 'generated-data root');
  const index = reports.map(({ channel, originalSlug, routeSlug, title, summary }) => ({
    channel,
    originalSlug,
    routeSlug,
    title,
    summary,
  }));
  await writeFileAtomically(outputRoot, generatedDataFile, JSON.stringify(index));
}

async function writeGeneratedReportFiles(reportsRoot, reports) {
  const safeReportsRoot = await ensureRealDirectory(reportsRoot, 'reports root');
  for (const report of reports) {
    const channel = safePathSegment(report.channel, 'report channel').toLowerCase();
    const routeSlug = safePathSegment(report.routeSlug, 'report route slug');
    const channelDir = resolveWithin(safeReportsRoot, channel);
    await assertWritableWithin(safeReportsRoot, channelDir, 'report channel directory');
    await ensureDir(channelDir);
    await writeFileAtomically(
      safeReportsRoot,
      resolveWithin(channelDir, `${routeSlug}.json`),
      JSON.stringify(report),
    );
  }
}

async function cleanGeneratedData(generatedDataFile, reportsRoot) {
  const safeDataFile = await assertRemovableWithin(GENERATED_DATA_ROOT, generatedDataFile, 'generated data file');
  const safeReportsRoot = await assertRemovableWithin(GENERATED_DATA_ROOT, reportsRoot, 'reports root');
  await fs.rm(safeDataFile, { force: true });
  await fs.rm(safeReportsRoot, { force: true, recursive: true });
}

function validateGeneratorDataset(data) {
  if (!data || typeof data !== 'object' || Array.isArray(data)) {
    throw new Error('Supermodel data must be an object keyed by channel');
  }

  const channelKeys = Object.keys(data).sort((a, b) => a.localeCompare(b));
  if (channelKeys.length === 0) {
    throw new Error('Supermodel data must contain at least one channel');
  }
  for (const channel of channelKeys) {
    safePathSegment(channel, 'channel');
    const models = data[channel];
    if (!models || typeof models !== 'object' || Array.isArray(models)) {
      throw new Error(`Channel ${channel} must contain an object keyed by model`);
    }
    const modelEntries = Object.entries(models);
    if (modelEntries.length === 0) {
      throw new Error(`Channel ${channel} must contain at least one model`);
    }
    for (const [modelName, modelData] of modelEntries) {
      safePathSegment(modelName, 'model name');
      if (!modelData || typeof modelData !== 'object' || Array.isArray(modelData)) {
        throw new Error(`Model ${channel}/${modelName} must contain an object`);
      }
    }
  }
  return channelKeys;
}

async function main() {
  const options = parseArgs();
  const data = runPythonExtraction(options.source);
  const channelKeys = validateGeneratorDataset(data);

  if (options.cleanData) {
    await cleanGeneratedData(options.generatedDataFile, options.reportsRoot);
  }

  options.markdownDir = await ensureRealDirectory(options.markdownDir, 'markdown root');

  const slugTracker = new Map();
  const reports = [];
  let markdownCount = 0;

  const totalTasks = channelKeys.reduce((sum, channel) => {
    const models = data[channel];
    if (!models || typeof models !== 'object') return sum;
    return sum + Object.keys(models).length;
  }, 0);
  const progress = createProgressTracker(totalTasks, 'Generating omnimodel markdown');
  let processedTasks = 0;

  for (const channel of channelKeys) {
    const safeChannel = safePathSegment(channel, 'channel');
    const models = data[channel];
    const modelEntries = Object.entries(models).sort((a, b) => a[0].localeCompare(b[0]));
    for (const [modelName, modelData] of modelEntries) {
      const safeModelName = safePathSegment(modelName, 'model name');
      const routeSlugBase = sanitizeRouteSlug(modelName);
      const routeSlug = ensureUniqueSlug(slugTracker, channel, routeSlugBase);
      const report = buildReport(channel, modelName, routeSlug, modelData);
      reports.push(report);

      const markdownDir = resolveWithin(options.markdownDir, safeChannel, safeModelName);
      await assertWritableWithin(options.markdownDir, markdownDir, 'markdown model directory');
      const imageSection = await collectImageMarkdown(markdownDir);
      const markdownContent = buildMarkdown(channel, modelName, modelData, imageSection);
      await writeMarkdownFile(options.markdownDir, channel, modelName, markdownContent);
      markdownCount += 1;
      processedTasks += 1;
      progress.update(processedTasks);
    }
  }

  reports.sort((a, b) => {
    const channelCompare = a.channel.localeCompare(b.channel);
    if (channelCompare !== 0) return channelCompare;
    return a.originalSlug.localeCompare(b.originalSlug);
  });

  await writeGeneratedDataFile(options.generatedDataFile, reports);
  await writeGeneratedReportFiles(options.reportsRoot, reports);

  if (totalTasks === 0) {
    console.log('No omnimodel markdown to generate.');
  } else {
    progress.update(totalTasks);
    console.log(`Generated ${markdownCount} markdown files.`);
  }
  console.log(`Generated data file: ${options.generatedDataFile}`);
  console.log(`Generated report files: ${options.reportsRoot}`);
}

main().catch((error) => {
  console.error('Failed to generate omnimodel outputs:');
  console.error(error);
  process.exit(1);
});

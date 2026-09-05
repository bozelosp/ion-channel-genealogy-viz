#!/usr/bin/env node

import { spawn, spawnSync } from 'node:child_process';
import { promises as fs } from 'node:fs';
import path from 'node:path';
import { fileURLToPath } from 'node:url';
import {
  assertWritableWithin,
  ensureRealDirectory,
  resolveWithin,
  safePathSegment,
} from './path-safety.mjs';

const __filename = fileURLToPath(import.meta.url);
const __dirname = path.dirname(__filename);
const projectRoot = path.resolve(__dirname, '..');
const PICKLE_LOADER = path.join(__dirname, 'load-supermodel-json.py');

function printUsage() {
  console.log(`Usage: node scripts/render-omnimodel-plots.mjs [options]

Options:
  --source <path>         Path to supermodel_data.pkl (default: supermodels/supermodel_data.pkl)
  --out-markdown <path>   Markdown root (default: public/data/omnimodels)
  --states <list>         Comma-separated list of gates to render (e.g., m,h)
  --overwrite             Re-render plots even if files exist
  --dry-run               Preview which plots would be generated
  --help                  Show this help message
`);
}

function parseArgs() {
  const args = process.argv.slice(2);
  const options = {
    source: path.resolve(projectRoot, 'supermodels', 'supermodel_data.pkl'),
    markdownDir: path.join(projectRoot, 'public', 'data', 'omnimodels'),
    statesFilter: null,
    overwrite: false,
    dryRun: false,
  };

  for (let i = 0; i < args.length; i += 1) {
    const arg = args[i];
    if (arg === '--source' && args[i + 1]) {
      options.source = path.resolve(args[i + 1]);
      i += 1;
    } else if (arg === '--out-markdown' && args[i + 1]) {
      options.markdownDir = path.resolve(args[i + 1]);
      i += 1;
    } else if (arg === '--states' && args[i + 1]) {
      options.statesFilter = new Set(
        args[i + 1]
          .split(',')
          .map((token) => token.trim())
          .filter(Boolean),
      );
      i += 1;
    } else if (arg === '--overwrite') {
      options.overwrite = true;
    } else if (arg === '--dry-run') {
      options.dryRun = true;
    } else if (arg === '--help' || arg === '-h') {
      printUsage();
      process.exit(0);
    } else {
      throw new Error(`Unknown argument: ${arg}`);
    }
  }

  return options;
}

function pathExists(p) {
  return fs
    .access(p)
    .then(() => true)
    .catch(() => false);
}

function sanitizeLabel(rawLabel) {
  if (!rawLabel) return 'baseline';
  const trimmed = rawLabel.trim();
  if (!trimmed) return 'baseline';
  const numeric = Number(trimmed);
  if (!Number.isNaN(numeric)) {
    return `${numeric}\u00b0C`;
  }
  return trimmed.replace(/_/g, ' ');
}

function splitSeriesByState(seriesMap) {
  const grouped = new Map();
  if (!seriesMap || typeof seriesMap !== 'object') return grouped;
  for (const [key, values] of Object.entries(seriesMap)) {
    if (!Array.isArray(values) || values.length === 0) continue;
    const [base, suffix] = key.split('_', 2);
    if (!grouped.has(base)) {
      grouped.set(base, []);
    }
    grouped.get(base).push({ label: sanitizeLabel(suffix), values });
  }
  return grouped;
}

function formatTitle(channel, modelName, stateName) {
  return `${modelName} – ${channel} – ${stateName} gate`;
}

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

function runPythonRenderer(payload) {
  const result = spawnSync('python3', [PICKLE_LOADER, payload.source], {
    encoding: 'utf-8',
    maxBuffer: 128 * 1024 * 1024,
    timeout: 60_000,
  });

  if (result.error) {
    throw result.error;
  }

  if (result.status !== 0) {
    throw new Error(`Failed to load pickle: ${result.stderr}`);
  }

  return JSON.parse(result.stdout);
}

function runPythonPlotter(payload) {
  const pythonCode = `import json, math, os, sys, tempfile, time
import pathlib

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.ticker import ScalarFormatter

THEME_FACE = '#f8fafc'
THEME_GRID = {
    'color': '#cbd5f5',
    'linewidth': 0.7,
    'linestyle': '--',
    'alpha': 0.6,
}

plt.rcParams.update({
    'figure.dpi': 120,
    'savefig.dpi': 180,
    'axes.edgecolor': '#94a3b8',
    'axes.labelcolor': '#0f172a',
    'axes.titlesize': 14,
    'axes.labelsize': 12,
    'xtick.color': '#334155',
    'ytick.color': '#334155',
    'font.family': 'DejaVu Sans',
})

def format_seconds(value):
    if value is None or not math.isfinite(value):
        return '?:??'
    value = max(0, int(round(value)))
    minutes, seconds = divmod(value, 60)
    if minutes >= 60:
        hours, minutes = divmod(minutes, 60)
        return f"{hours}h {minutes}m"
    return f"{minutes}:{seconds:02d}"

def render_state(fig_cfg):
    voltage = fig_cfg['voltage']
    steady_series = fig_cfg['steady']
    tau_series = fig_cfg['tau']
    outfile = pathlib.Path(fig_cfg['output'])
    outfile.parent.mkdir(parents=True, exist_ok=True)

    fig, (ax_ss, ax_tau) = plt.subplots(1, 2, figsize=fig_cfg['figsize'], constrained_layout=True)
    for ax in (ax_ss, ax_tau):
        ax.set_facecolor(THEME_FACE)
        ax.grid(**THEME_GRID)
        for spine in ax.spines.values():
            spine.set_alpha(0.4)

    for series in steady_series:
        ax_ss.plot(
            voltage,
            series['values'],
            label=series['label'],
            linewidth=2.3,
            solid_capstyle='round',
        )

    ax_ss.set_title('Steady-state activation')
    ax_ss.set_xlabel('Membrane voltage (mV)')
    ax_ss.set_ylabel('Probability')
    ax_ss.set_ylim(-0.05, 1.05)

    tau_y_min = math.inf
    tau_y_max = 0
    for series in tau_series:
        values = series['values']
        ax_tau.plot(
            voltage,
            values,
            label=series['label'],
            linewidth=2.3,
            solid_capstyle='round',
        )
        tau_y_min = min(tau_y_min, min(values))
        tau_y_max = max(tau_y_max, max(values))

    ax_tau.set_title('Time constant (ms)')
    ax_tau.set_xlabel('Membrane voltage (mV)')
    ax_tau.set_ylabel('Tau (ms)')

    if tau_y_max / max(tau_y_min, 1e-6) > 25:
        ax_tau.set_yscale('log')
        ax_tau.yaxis.set_major_formatter(ScalarFormatter())
    else:
        ax_tau.set_ylim(bottom=0)

    handles, labels = ax_ss.get_legend_handles_labels()
    if not handles:
        handles, labels = ax_tau.get_legend_handles_labels()
    if handles:
        fig.legend(
            handles,
            labels,
            loc='upper center',
            bbox_to_anchor=(0.5, 1.02),
            ncol=min(3, len(labels)),
            frameon=False,
            fontsize=10,
        )

    fig.suptitle(fig_cfg['title'], fontsize=16, color='#0f172a', y=1.05)
    descriptor, temporary = tempfile.mkstemp(
        prefix=f'.{outfile.name}.',
        suffix='.tmp.png',
        dir=outfile.parent,
    )
    os.close(descriptor)
    try:
        fig.savefig(temporary, format='png', bbox_inches='tight', facecolor='white')
        os.chmod(temporary, 0o644)
        os.replace(temporary, outfile)
    finally:
        try:
            os.unlink(temporary)
        except FileNotFoundError:
            pass
        plt.close(fig)

config = json.loads(sys.stdin.read())
states = config['states']
start = time.time()
for idx, state_cfg in enumerate(states, start=1):
    render_state(state_cfg)
    elapsed = max(time.time() - start, 1e-6)
    rate = idx / elapsed
    remaining = len(states) - idx
    eta = remaining / rate if rate > 0 else 0
    sys.stdout.write(f"Rendering omnimodel plots: {idx}/{len(states)} (ETA {format_seconds(eta)})\r")
    sys.stdout.flush()
if states:
    sys.stdout.write("\n")
`;

  return new Promise((resolve, reject) => {
    const child = spawn('python3', ['-c', pythonCode], {
      stdio: ['pipe', 'pipe', 'inherit'],
    });

    child.on('error', reject);
    child.stdout.on('data', (data) => {
      process.stdout.write(data);
    });
    child.on('close', (code) => {
      if (code === 0) {
        resolve();
      } else {
        reject(new Error(`Python renderer failed with exit code ${code}`));
      }
    });

    child.stdin.write(JSON.stringify(payload));
    child.stdin.end();
  });
}

async function main() {
  const options = parseArgs();
  const dataset = runPythonRenderer({ source: options.source });
  const markdownRoot = await ensureRealDirectory(options.markdownDir, 'markdown root');

  const targetStates = options.statesFilter;
  const plannedStates = [];
  const pythonPayload = { states: [] };

  for (const [channel, models] of Object.entries(dataset)) {
    if (!models || typeof models !== 'object') continue;
    const safeChannel = safePathSegment(channel, 'channel');

    for (const [modelName, modelData] of Object.entries(models)) {
      const safeModelName = safePathSegment(modelName, 'model name');
      const voltageSeries = modelData?.RATE_VALS_V ?? {};
      const steadySeries = splitSeriesByState(modelData?.RATE_VALS_SS ?? {});
      const tauSeries = splitSeriesByState(modelData?.RATE_VALS_TAU ?? {});

      for (const [stateName, voltage] of Object.entries(voltageSeries)) {
        if (!Array.isArray(voltage) || voltage.length === 0) continue;
        if (targetStates && !targetStates.has(stateName)) continue;
        const safeStateName = safePathSegment(stateName, 'state name');

        const steady = steadySeries.get(stateName) ?? [];
        const tau = tauSeries.get(stateName) ?? [];
        if (steady.length === 0 && tau.length === 0) continue;

        const markdownDir = resolveWithin(markdownRoot, safeChannel, safeModelName);
        const imagesDir = resolveWithin(markdownDir, 'images');
        const outputPath = resolveWithin(imagesDir, `${safeStateName}_summary.png`);
        await assertWritableWithin(markdownRoot, imagesDir, 'plot output directory');
        if (!options.dryRun) {
          await fs.mkdir(imagesDir, { recursive: true });
        }
        await assertWritableWithin(markdownRoot, outputPath, 'plot output file');

        const destinationExists = await pathExists(outputPath);
        if (destinationExists && !options.overwrite) {
          continue;
        }

        plannedStates.push({ channel, modelName, stateName, outputPath });
        pythonPayload.states.push({
          title: formatTitle(channel, modelName, stateName),
          output: outputPath,
          figsize: [8.5, 4.5],
          voltage,
          steady,
          tau,
        });
      }
    }
  }

  if (options.dryRun) {
    if (plannedStates.length === 0) {
      console.log('No plots queued (dry run).');
    } else {
      console.log('Planned omnimodel plots:');
      plannedStates.forEach((item) => {
        console.log(` - ${item.channel}/${item.modelName} [${item.stateName}] -> ${item.outputPath}`);
      });
    }
    return;
  }

  if (pythonPayload.states.length === 0) {
    console.log('No plots queued (existing files up to date or no matching states).');
    return;
  }

  const queueProgress = createProgressTracker(plannedStates.length, 'Queuing omnimodel plots');
  plannedStates.forEach((_, index) => queueProgress.update(index + 1));
  queueProgress.update(plannedStates.length);

  await runPythonPlotter(pythonPayload);
  console.log(`Rendered ${pythonPayload.states.length} omnimodel plot(s).`);
}

main().catch((error) => {
  console.error('Failed to render omnimodel plots:');
  console.error(error.message ?? error);
  process.exit(1);
});

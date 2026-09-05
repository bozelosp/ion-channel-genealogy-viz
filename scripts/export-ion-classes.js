#!/usr/bin/env node
/*
  Export per-ion-class summaries from public/data/network_data.json

  Outputs:
    - exports/ion-classes/csv/<IonClass>.csv
    - exports/ion-classes/json/<IonClass>.json

  CSV columns:
    node_id,name,ion_class,copies,unique_ids
      - copies = 1 + (identical_models.length or num_of_identicals if > 0)
      - unique_ids = semicolon-separated list including the node's unique ID and all identicals' unique IDs

  JSON format per file:
    {
      ion_class: "K" | "Na" | "Ca" | "Ih" | "KCa" | "Other",
      nodes_count: number,
      copies_total: number,
      nodes: [
        { node_id, name, ion_class, copies, unique_ids: string[] }
      ]
    }
*/

const fs = require('fs');
const path = require('path');

const ROOT = process.cwd();
const INPUT = path.join(ROOT, 'public', 'data', 'network_data.json');
const OUT_BASE = path.join(ROOT, 'exports', 'ion-classes');
const OUT_CSV = path.join(OUT_BASE, 'csv');
const OUT_JSON = path.join(OUT_BASE, 'json');

const CANONICAL_CLASSES = ['K', 'Na', 'Ca', 'Ih', 'KCa', 'Other'];

function canonicalClass(raw) {
  if (!raw) return 'Other';
  const s = String(raw).trim();
  if (!s) return 'Other';
  if (s.toLowerCase() === 'ih' || s.toLowerCase() === 'h') return 'Ih';
  if (s.toLowerCase() === 'kca') return 'KCa';
  const upper = s[0].toUpperCase() + s.slice(1);
  if (CANONICAL_CLASSES.includes(upper)) return upper;
  return 'Other';
}

function ensureDir(dir) {
  fs.mkdirSync(dir, { recursive: true });
}

function escCSV(v) {
  const s = (v ?? '').toString();
  return /[",\n]/.test(s) ? '"' + s.replace(/"/g, '""') + '"' : s;
}

function main() {
  if (!fs.existsSync(INPUT)) {
    console.error(`Input not found: ${INPUT}`);
    process.exit(1);
  }
  const raw = fs.readFileSync(INPUT, 'utf8');
  let data;
  try {
    data = JSON.parse(raw);
  } catch (e) {
    console.error('Failed to parse JSON:', e.message);
    process.exit(1);
  }

  const nodes = Array.isArray(data?.nodes) ? data.nodes : [];

  // Group nodes by ion class
  const byClass = new Map();
  for (const node of nodes) {
    // Determine class
    const ionClass = canonicalClass(node?.original_model?.ion_class || node?.ion_class);
    if (!byClass.has(ionClass)) byClass.set(ionClass, []);

    // Extract id/name
    const node_id = String(node?.id ?? '');
    const name = String(node?.name ?? node?.label ?? node_id);

    // copies = 1 + identicals.length (preferred) or 1 + max(num_of_identicals, 0)
    const identList = Array.isArray(node?.identical_models) ? node.identical_models : [];
    const nFromField = typeof node?.num_of_identicals === 'number' ? Math.max(0, node.num_of_identicals) : 0;
    const copies = 1 + (identList.length || nFromField);

    // unique ids: include primary and all identicals
    const ids = [];
    const primary = node?.original_model?.unique_modelDB_mod_id;
    if (primary) ids.push(String(primary));
    for (const e of identList) {
      const id = e?.unique_modelDB_mod_id;
      if (id) ids.push(String(id));
    }
    // Dedupe while preserving order
    const seen = new Set();
    const unique_ids = ids.filter((x) => (x && !seen.has(x) && seen.add(x)));

    byClass.get(ionClass).push({ node_id, name, ion_class: ionClass, copies, unique_ids });
  }

  ensureDir(OUT_CSV);
  ensureDir(OUT_JSON);

  const summary = [];

  for (const ionClass of CANONICAL_CLASSES) {
    const rows = byClass.get(ionClass) || [];
    const nodes_count = rows.length;
    const copies_total = rows.reduce((acc, r) => acc + (r.copies || 0), 0);
    summary.push({ ion_class: ionClass, nodes_count, copies_total });

    // CSV
    const csvPath = path.join(OUT_CSV, `${ionClass}.csv`);
    const csvLines = [];
    csvLines.push(`# ion_class=${ionClass}, nodes=${nodes_count}, copies_total=${copies_total}`);
    csvLines.push(['node_id','name','ion_class','copies','unique_ids'].map(escCSV).join(','));
    for (const r of rows) {
      const uniqueJoined = r.unique_ids.join(';');
      csvLines.push([r.node_id, r.name, r.ion_class, String(r.copies), uniqueJoined].map(escCSV).join(','));
    }
    fs.writeFileSync(csvPath, csvLines.join('\n'));

    // JSON
    const jsonPath = path.join(OUT_JSON, `${ionClass}.json`);
    const jsonObj = { ion_class: ionClass, nodes_count, copies_total, nodes: rows };
    fs.writeFileSync(jsonPath, JSON.stringify(jsonObj, null, 2));
  }

  // Write overall summary
  const sumPath = path.join(OUT_BASE, 'SUMMARY.json');
  fs.writeFileSync(sumPath, JSON.stringify({ generated_at: new Date().toISOString(), summary }, null, 2));

  console.log('Export complete.');
  for (const s of summary) {
    console.log(`${s.ion_class}: nodes=${s.nodes_count}, copies_total=${s.copies_total}`);
  }
  console.log(`\nCSV:  ${OUT_CSV}`);
  console.log(`JSON: ${OUT_JSON}`);
}

main();


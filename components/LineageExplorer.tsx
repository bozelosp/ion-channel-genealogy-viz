'use client';

import Link from 'next/link';
import { useRef, useState, type KeyboardEvent } from 'react';
import { comparisonPath } from '@/lib/comparison';
import type { Lineage, LineageLayout, LineageNode } from '@/lib/lineage';

// Node colours follow the graph view: gold origins, sky-blue models.
const DOT = {
  origin: 'bg-[#ffd700] ring-1 ring-[#215885] dark:bg-[#fbbf24] dark:ring-[#1e40af]',
  close: 'bg-[#00BFFF]/35 ring-1 ring-[#00BFFF] dark:bg-[#3b82f6]/40 dark:ring-[#3b82f6]',
  edited: 'bg-[#00BFFF] dark:bg-[#3b82f6]',
  distant: 'bg-[#215885] dark:bg-blue-200',
};

function dotClass(node: LineageNode): string {
  if (node.parent < 0) return DOT.origin;
  const sim = node.simParent ?? 0;
  return sim >= 99 ? DOT.close : sim >= 95 ? DOT.edited : DOT.distant;
}

const percent = (value: number, digits: number) => value.toFixed(digits) + '%';
const formatThreshold = (value: number) => (Number.isInteger(value) ? String(value) : value.toFixed(2));

export default function LineageExplorer({
  lineage,
  layout,
  initialId,
}: {
  lineage: Lineage;
  layout: LineageLayout;
  initialId: string;
}) {
  const { nodes, file } = lineage;
  const linkThreshold = formatThreshold(lineage.linkThreshold);
  const [selected, setSelected] = useState(() => Math.max(0, nodes.findIndex((node) => node.id === initialId)));
  const buttons = useRef<(HTMLButtonElement | null)[]>([]);

  const chain: number[] = [];
  for (let i = selected; i >= 0; i = nodes[i].parent) chain.unshift(i);
  const inChain = new Set(chain.slice(1));
  const restEdges = layout.edges.filter((_, i) => !inChain.has(i)).join('');
  const chainEdges = chain.slice(1).map((i) => layout.edges[i]).join('');

  const labelled = new Set(nodes.flatMap((node, i) => (node.parent < 0 || node.copies >= 3 ? [i] : [])));
  labelled.add(selected);

  const node = nodes[selected];
  const parent = node.parent >= 0 ? nodes[node.parent] : null;
  const origin = nodes[node.origin];
  const graphQuery = new URLSearchParams({ q: node.id });
  if (node.copies < 2) graphQuery.set('copies', '1');

  const totalCopies = nodes.reduce((sum, n) => sum + n.copies, 0);
  const origins = nodes.filter((n) => n.parent < 0).length;
  const toPercent = (value: number) => (value / layout.size) * 100 + '%';

  function select(index: number) {
    setSelected(index);
    buttons.current[index]?.focus();
  }

  function onKeyDown(event: KeyboardEvent<HTMLButtonElement>) {
    const moves: Record<string, number> = {
      ArrowRight: selected + 1,
      ArrowDown: selected + 1,
      ArrowLeft: selected - 1,
      ArrowUp: selected - 1,
      Home: 0,
      End: nodes.length - 1,
    };
    const target = moves[event.key];
    if (target === undefined) return;
    event.preventDefault();
    select((target + nodes.length) % nodes.length);
  }

  return (
    <figure className="m-0 overflow-hidden rounded-xl border border-slate-200 bg-white dark:border-slate-800 dark:bg-slate-900">
      <figcaption className="flex flex-wrap items-baseline justify-between gap-x-4 gap-y-1 border-b border-slate-200 px-4 py-3 dark:border-slate-800">
        <span className="flex flex-wrap items-baseline gap-x-3 gap-y-1">
          <span className="font-mono text-sm font-medium">{file}</span>
          <span className="text-sm text-slate-600 dark:text-slate-400">
            {nodes.length} variants · {totalCopies} identical copies · {origins} independent origins
          </span>
        </span>
        <span className="text-sm text-slate-600 dark:text-slate-400">Radius: year of ModelDB entry</span>
      </figcaption>

      <div className="flex flex-wrap">
        <div className="flex min-w-0 flex-[1_1_540px] justify-center border-b border-slate-200 p-3 dark:border-slate-800 lg:border-b-0 lg:border-r">
          <div role="radiogroup" aria-label={`Variants of ${file}`} className="relative aspect-square w-full max-w-[540px]">
            <svg viewBox={`0 0 ${layout.size} ${layout.size}`} aria-hidden="true" className="absolute inset-0 h-full w-full">
              {layout.rings.map((ring) => (
                <g key={ring.year}>
                  <circle cx={layout.size / 2} cy={layout.size / 2} r={ring.r} className="fill-none stroke-slate-200 dark:stroke-slate-800" />
                  <text x={layout.size / 2} y={layout.size / 2 - ring.r - 4} textAnchor="middle" className="fill-slate-400 font-mono text-[10.5px] dark:fill-slate-500">
                    {ring.year}
                  </text>
                </g>
              ))}
              <path d={restEdges} className="fill-none stroke-slate-400 dark:stroke-slate-500" strokeWidth={1.2} />
              <path d={chainEdges} className="fill-none stroke-blue-600 dark:stroke-blue-400" strokeWidth={2.2} strokeLinecap="round" />
            </svg>
            {nodes.map((n, i) => {
              const { x, y } = layout.points[i];
              const size = 2 * (2.5 + 1.25 * Math.sqrt(n.copies));
              const isSelected = i === selected;
              return (
                <button
                  key={n.id}
                  ref={(element) => { buttons.current[i] = element; }}
                  type="button"
                  role="radio"
                  aria-checked={isSelected}
                  aria-label={`${n.id}, ${n.year}`}
                  tabIndex={isSelected ? 0 : -1}
                  onClick={() => select(i)}
                  onKeyDown={onKeyDown}
                  className="absolute flex -translate-x-1/2 -translate-y-1/2 cursor-pointer items-center justify-center rounded-full focus-visible:outline-2 focus-visible:outline-offset-1 focus-visible:outline-blue-600 dark:focus-visible:outline-blue-400"
                  style={{ left: toPercent(x), top: toPercent(y), width: Math.max(24, size + 6), height: Math.max(24, size + 6) }}
                >
                  <span
                    className={`block shrink-0 rounded-full ${dotClass(n)} ${isSelected ? 'outline-2 outline-offset-2 outline-blue-600 dark:outline-blue-400' : ''}`}
                    style={{ width: size, height: size }}
                  />
                </button>
              );
            })}
            {[...labelled].map((i) => {
              const n = nodes[i];
              const { x, y } = layout.points[i];
              const centre = layout.size / 2;
              // Labels sit across the radial direction, clear of the edges to
              // children further out: above nodes near the horizontal axis,
              // beside nodes near the vertical axis.
              // Offsets are in CSS pixels from the dot's edge, so labels keep
              // their spacing when the figure scales down.
              const angle = Math.atan2(y - centre, x - centre);
              const above = i !== 0 && Math.abs(Math.cos(angle)) > 0.7;
              const left = i !== 0 && !above && x < centre;
              const gap = 2.5 + 1.25 * Math.sqrt(n.copies) + (i === selected ? 7 : 4);
              const transform = i === 0
                ? `translate(${gap * 0.7}px, ${gap * 0.7}px)`
                : above
                  ? `translate(-50%, calc(-100% - ${gap}px))`
                  : left
                    ? `translate(calc(-100% - ${gap}px), -50%)`
                    : `translate(${gap}px, -50%)`;
              // Phones keep only the origin and selected labels.
              const visibility = i === 0 || i === selected ? '' : 'max-sm:hidden';
              return (
                <span
                  key={n.id}
                  aria-hidden="true"
                  className={`pointer-events-none absolute whitespace-nowrap font-mono text-[11px] leading-[14px] ${visibility} ${i === selected ? 'text-blue-700 dark:text-blue-300' : 'text-slate-600 dark:text-slate-300'}`}
                  style={{ left: toPercent(x), top: toPercent(y), transform }}
                >
                  {n.name} · {n.year}
                </span>
              );
            })}
          </div>
        </div>

        <div className="flex min-w-0 flex-[1_1_340px] flex-col gap-5 p-5">
          <div className="flex flex-col gap-0.5">
            <span className="break-all font-mono text-xl font-medium">{node.id}</span>
            <span className="text-sm text-slate-600 dark:text-slate-400">
              <a href={`https://modeldb.science/${node.modelDb}`} className="underline-offset-2 hover:underline">ModelDB {node.modelDb}</a> · {node.year}
            </span>
          </div>

          <dl className="m-0 grid grid-cols-2 gap-px overflow-hidden rounded-lg border border-slate-200 bg-slate-200 dark:border-slate-800 dark:bg-slate-800">
            {[
              ['Similarity to parent', parent && node.simParent !== null ? percent(node.simParent, 2) : 'origin'],
              [`Similarity to ${origin.name}`, node.parent < 0 ? '—' : node.simOrigin === null ? `< ${linkThreshold}%` : percent(node.simOrigin, 1)],
              ['Identical copies', String(node.copies)],
              ['Modified descendants', String(node.descendants)],
            ].map(([term, value]) => (
              <div key={term} className="bg-white px-3.5 py-3 dark:bg-slate-900">
                <dt className="text-xs text-slate-600 dark:text-slate-400">{term}</dt>
                <dd className="m-0 mt-0.5 font-mono text-base">{value}</dd>
              </div>
            ))}
          </dl>

          <div className="flex flex-col gap-2">
            <span className="text-xs text-slate-600 dark:text-slate-400">Path from origin</span>
            <ol className="m-0 flex list-none flex-wrap items-center gap-1.5 p-0 font-mono text-xs">
              {chain.map((index, step) => (
                <li key={index} className="flex items-center gap-1.5">
                  {step > 0 && <span className="text-slate-500 dark:text-slate-400">→ {percent(nodes[index].simParent ?? 0, 2)} →</span>}
                  <span className={`rounded-md px-2 py-0.5 ${index === selected ? 'bg-blue-600 text-white' : 'bg-slate-100 text-slate-800 dark:bg-slate-800 dark:text-slate-100'}`}>
                    {nodes[index].name}
                  </span>
                </li>
              ))}
            </ol>
          </div>

          <div className="mt-auto flex flex-wrap gap-2.5">
            {parent && (
              <Link href={comparisonPath(parent.id, node.id)} className="inline-flex h-11 items-center rounded-lg bg-blue-600 px-4 text-sm font-semibold text-white transition-colors hover:bg-blue-700">
                Compare with parent
              </Link>
            )}
            <Link href={`/visualizer?${graphQuery}`} className="inline-flex h-11 items-center rounded-lg border border-slate-300 px-4 text-sm font-medium transition-colors hover:bg-slate-100 dark:border-slate-700 dark:hover:bg-slate-800">
              Show in graph
            </Link>
          </div>

          <ul className="m-0 flex list-none flex-wrap gap-x-3.5 gap-y-1.5 border-t border-slate-200 p-0 pt-3.5 text-xs text-slate-600 dark:border-slate-800 dark:text-slate-400">
            <li className="flex items-center gap-1.5"><span className={`h-2.5 w-2.5 rounded-full ${DOT.origin}`} />origin</li>
            <li className="flex items-center gap-1.5"><span className={`h-2.5 w-2.5 rounded-full ${DOT.close}`} />≥ 99% to parent</li>
            <li className="flex items-center gap-1.5"><span className={`h-2.5 w-2.5 rounded-full ${DOT.edited}`} />95–99%</li>
            <li className="flex items-center gap-1.5"><span className={`h-2.5 w-2.5 rounded-full ${DOT.distant}`} />&lt; 95%</li>
            <li>size: identical copies</li>
          </ul>
        </div>
      </div>

      <p className="m-0 border-t border-slate-200 px-4 py-2.5 text-xs text-slate-600 dark:border-slate-800 dark:text-slate-400">
        Each file is linked to its most similar earlier file (source similarity ≥ {linkThreshold}%). Undated and unlinked files are omitted.
      </p>
    </figure>
  );
}

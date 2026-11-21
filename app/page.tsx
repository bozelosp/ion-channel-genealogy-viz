import Link from 'next/link';
import { readFile } from 'node:fs/promises';
import path from 'node:path';
import omnimodelIndex from '@/public/data/omnimodel-index.json';
import LineageExplorer from '@/components/LineageExplorer';
import { comparisonPath } from '@/lib/comparison';
import { getVisibleIonClassOptions } from '@/lib/ionClasses';
import { buildLineage, layoutLineage } from '@/lib/lineage';

const LINEAGE_FILE = 'kaprox.mod';
const LINEAGE_ORIGIN = '2796_kaprox-ID-1';
const LINEAGE_SELECTED = '101629_kaprox-ID-1';
const EXAMPLE_COMPARISON = comparisonPath('2796_kaprox-ID-1', '101629_kaprox-ID-1');
const EXAMPLE_QUERIES = ['kaprox.mod', '149100', 'h.mod', 'na.mod'];

// Read once per process at build time; the page is prerendered.
let networkData: Promise<unknown> | undefined;
function loadNetworkData(): Promise<unknown> {
  networkData ??= readFile(path.join(process.cwd(), 'public', 'data', 'network_data.json'), 'utf8').then(JSON.parse);
  return networkData;
}

function count(value: unknown, label: string): number {
  if (!Array.isArray(value) || value.length === 0) throw new Error(`Missing ${label} in network data`);
  return value.length;
}

const format = (value: number) => value.toLocaleString('en-US');

export default async function Home() {
  const network = (await loadNetworkData()) as { nodes?: unknown; links?: unknown };
  const lineage = buildLineage(network, LINEAGE_FILE, LINEAGE_ORIGIN);
  const layout = layoutLineage(lineage);
  const fileCount = count(network.nodes, 'nodes');
  const linkCount = count(network.links, 'links');
  const reportChannels = new Set(omnimodelIndex.map((report) => report.channel.toLowerCase()));
  const channels = getVisibleIonClassOptions().filter((channel) => reportChannels.has(channel.toLowerCase()));

  return (
    <div className="min-h-screen bg-slate-50 font-sans text-slate-900 dark:bg-slate-950 dark:text-slate-100">
      <div className="mx-auto max-w-[1200px] px-4 sm:px-6 lg:px-8">
        <header className="flex flex-wrap items-center justify-between gap-x-6 gap-y-3 border-b border-slate-200 py-4 dark:border-slate-800">
          <Link href="/" className="flex items-center gap-2.5 text-[17px] font-semibold tracking-tight">
            <svg width="24" height="24" viewBox="0 0 24 24" aria-hidden="true" className="fill-none stroke-[#00BFFF] dark:stroke-[#3b82f6]" strokeWidth={1.6} strokeLinecap="round">
              <path d="M5 12h4M9 12l5-6M9 12l5 6M17 6.5l2.6 4.2" />
              <circle cx="5" cy="12" r="2.4" className="fill-[#ffd700] stroke-none dark:fill-[#fbbf24]" />
              <circle cx="15.5" cy="5" r="2" className="fill-[#00BFFF] stroke-none dark:fill-[#3b82f6]" />
              <circle cx="15.5" cy="19" r="2" className="fill-[#00BFFF] stroke-none dark:fill-[#3b82f6]" />
              <circle cx="20.5" cy="12" r="1.6" className="fill-slate-400 stroke-none dark:fill-slate-500" />
            </svg>
            ICG Explorer
          </Link>
          <nav aria-label="Main" className="flex flex-wrap items-center gap-x-6 gap-y-1 text-sm text-slate-600 dark:text-slate-400">
            <Link href="/visualizer" className="hover:text-slate-900 dark:hover:text-slate-100">Graph</Link>
            <Link href="/omnimodels" className="hover:text-slate-900 dark:hover:text-slate-100">Kinetic reports</Link>
            <Link href={EXAMPLE_COMPARISON} className="hover:text-slate-900 dark:hover:text-slate-100">Compare</Link>
          </nav>
        </header>

        <main>
          <section className="flex flex-wrap items-end justify-between gap-x-12 gap-y-6 pb-7 pt-12">
            <div className="flex flex-[1_1_420px] flex-col gap-2.5">
              <h1 className="m-0 text-4xl font-semibold leading-tight tracking-tight sm:text-5xl">Ion Channel Genealogy</h1>
              <p className="m-0 max-w-[52ch] text-[17px] text-slate-600 dark:text-slate-400">
                Reuse, modification and kinetics of ion-channel models in ModelDB.
              </p>
            </div>
            <div className="flex max-w-[520px] flex-[1_1_420px] flex-col gap-2.5">
              <form action="/visualizer" method="get" role="search" className="flex h-[52px] items-center gap-2.5 rounded-xl border border-slate-300 bg-white pl-3.5 pr-1.5 focus-within:border-blue-500 focus-within:ring-2 focus-within:ring-blue-500 dark:border-slate-700 dark:bg-slate-900 dark:focus-within:border-blue-400 dark:focus-within:ring-blue-400">
                <svg width="16" height="16" viewBox="0 0 24 24" aria-hidden="true" className="shrink-0 fill-none stroke-slate-400 dark:stroke-slate-500" strokeWidth={2} strokeLinecap="round">
                  <circle cx="11" cy="11" r="7" />
                  <path d="M20 20l-4-4" />
                </svg>
                <label htmlFor="landing-search" className="sr-only">ModelDB ID or .mod file name</label>
                <input
                  id="landing-search"
                  name="q"
                  type="search"
                  required
                  minLength={3}
                  maxLength={128}
                  placeholder="ModelDB ID or .mod file"
                  className="min-w-0 flex-1 bg-transparent font-mono text-[15px] text-slate-900 placeholder-slate-400 outline-none dark:text-slate-100 dark:placeholder-slate-500"
                />
                <button type="submit" className="inline-flex h-10 shrink-0 items-center rounded-lg bg-blue-600 px-4 text-sm font-semibold text-white transition-colors hover:bg-blue-700">
                  Search
                </button>
              </form>
              <div className="flex flex-wrap items-center gap-x-2 gap-y-1.5 text-[13px] text-slate-600 dark:text-slate-400">
                <span>Examples</span>
                {EXAMPLE_QUERIES.map((query) => (
                  <Link
                    key={query}
                    href={`/visualizer?${new URLSearchParams({ q: query })}`}
                    className="rounded-md border border-slate-200 px-2 py-0.5 font-mono text-slate-700 hover:border-slate-400 dark:border-slate-800 dark:text-slate-300 dark:hover:border-slate-600"
                  >
                    {query}
                  </Link>
                ))}
              </div>
            </div>
          </section>

          <LineageExplorer lineage={lineage} layout={layout} initialId={LINEAGE_SELECTED} />

          <nav aria-label="Tools" className="mt-6 grid grid-cols-[repeat(auto-fit,minmax(min(300px,100%),1fr))] gap-px overflow-hidden rounded-xl border border-slate-200 bg-slate-200 dark:border-slate-800 dark:bg-slate-800">
            {[
              { href: '/visualizer', title: 'Genealogy graph', detail: `${format(fileCount)} .mod files · ${format(linkCount)} similarity links` },
              { href: '/omnimodels', title: 'Kinetic reports', detail: `${format(omnimodelIndex.length)} models · ${channels.join(', ')}` },
              { href: EXAMPLE_COMPARISON, title: 'Source comparison', detail: 'line-level diff of any two files' },
            ].map((tool) => (
              <Link key={tool.title} href={tool.href} className="flex flex-col gap-0.5 bg-slate-50 px-5 py-4 transition-colors hover:bg-white dark:bg-slate-950 dark:hover:bg-slate-900">
                <span className="text-base font-semibold">{tool.title} →</span>
                <span className="font-mono text-[12.5px] text-slate-600 dark:text-slate-400">{tool.detail}</span>
              </Link>
            ))}
          </nav>
        </main>

        <footer className="flex flex-wrap justify-between gap-x-8 gap-y-3 pb-10 pt-7 text-[13px] text-slate-600 dark:text-slate-400">
          <span>
            Data: <a href="https://modeldb.science" className="underline-offset-2 hover:underline">ModelDB</a> ·{' '}
            <a href="https://icg.neurotheory.ox.ac.uk/" className="underline-offset-2 hover:underline">Ion Channel Genealogy</a>
          </span>
        </footer>
      </div>
    </div>
  );
}

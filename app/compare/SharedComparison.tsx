'use client';

import Link from 'next/link';
import { useSearchParams } from 'next/navigation';
import { useEffect, useState } from 'react';
import { ComparisonLinkButton } from '@/components/ComparisonLinkButton';
import { fetchModelSource, isModelSourceId } from '@/lib/comparison';
import { generateDiffHtml } from '../visualizer/generateDiffHtml';

type SourcePair = [string, string];

function PairDiff({ left, right }: { left: string; right: string }) {
  const [files, setFiles] = useState<SourcePair | null>(null);
  const [error, setError] = useState('');
  const [html, setHtml] = useState('');
  const [full, setFull] = useState(false);

  useEffect(() => {
    const controller = new AbortController();
    Promise.all([fetchModelSource(left, controller.signal), fetchModelSource(right, controller.signal)])
      .then(setFiles)
      .catch((cause: Error) => {
        if (!controller.signal.aborted) setError(cause.message);
      });
    return () => controller.abort();
  }, [left, right]);

  useEffect(() => {
    if (!files) return;
    let active = true;
    void generateDiffHtml(files[0], files[1], 'light', !full).then((next) => {
      if (active) setHtml(next);
    });
    return () => { active = false; };
  }, [files, full]);

  return <>
    <div className="flex flex-wrap items-center justify-between gap-3 border-y border-stone-300 py-3">
      <label className="flex items-center gap-2 text-sm">
        <input type="checkbox" checked={full} onChange={(event) => setFull(event.target.checked)} />
        Show unchanged lines
      </label>
      {files && <ComparisonLinkButton left={left} right={right} />}
    </div>
    <div className="grid grid-cols-2 gap-3 py-4 font-mono text-xs sm:text-sm">
      {[left, right].map(id => <a className="break-all underline underline-offset-4" key={id} href={'/api/source-code/' + encodeURIComponent(id)}>{id}</a>)}
    </div>
    {error ? <p role="alert">{error}</p> : html
      ? <iframe title="Code diff" className="h-[65vh] min-h-80 w-full rounded border border-stone-300 bg-white" sandbox="" referrerPolicy="no-referrer" srcDoc={html} />
      : <p role="status">Loading comparison…</p>}
  </>;
}

export default function SharedComparison() {
  const params = useSearchParams();
  const left = params.get('left') ?? '';
  const right = params.get('right') ?? '';
  const valid = params.getAll('left').length === 1 && params.getAll('right').length === 1
    && isModelSourceId(left) && isModelSourceId(right) && left !== right;

  return <main className="min-h-screen bg-[#f8f7f2] px-4 py-6 text-[#182a25] sm:px-10">
    <div className="mx-auto max-w-6xl">
      <nav className="flex flex-wrap justify-between gap-4 text-sm">
        <Link href="/">Ion Channel Genealogy</Link><Link href="/visualizer">Explore the graph →</Link>
      </nav>
      <h1 className="my-6 font-serif text-4xl sm:text-5xl">Source comparison</h1>
      {valid ? <PairDiff key={left + '/' + right} left={left} right={right} />
        : <p role="alert">This comparison needs two distinct, valid model source IDs. Open the graph to choose a pair.</p>}
    </div>
  </main>;
}

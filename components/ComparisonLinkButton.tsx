'use client';

import { useState } from 'react';
import { comparisonPath } from '@/lib/comparison';

export function ComparisonLinkButton({ left, right }: { left?: string; right?: string }) {
  const [status, setStatus] = useState<'ready' | 'copied' | 'manual'>('ready');
  if (!left || !right) return null;
  const href = comparisonPath(left, right);
  const copy = async () => {
    try {
      await navigator.clipboard.writeText(new URL(href, window.location.origin).href);
      setStatus('copied');
    } catch {
      setStatus('manual');
    }
  };
  if (status === 'manual') {
    return <a className="underline" href={href}>Open comparison link to copy its URL</a>;
  }
  return (
    <button type="button" onClick={copy} className="rounded-md border border-current px-3 py-1.5 text-sm font-medium" aria-live="polite">
      {status === 'copied' ? 'Link copied' : 'Copy comparison link'}
    </button>
  );
}

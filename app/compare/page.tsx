import { Suspense } from 'react';
import type { Metadata } from 'next';
import SharedComparison from './SharedComparison';

export const metadata: Metadata = { title: 'Compare model source · Ion Channel Genealogy' };

export default function ComparePage() {
  return <Suspense fallback={<p className="p-6">Loading comparison…</p>}><SharedComparison /></Suspense>;
}

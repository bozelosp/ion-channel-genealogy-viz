"use client";

import Link from "next/link";
import { useDeferredValue, useEffect, useMemo, useState } from "react";

const PAGE_SIZE = 60;

interface ModelGridProps {
  channel: string;
  models: {
    slug: string;
    title: string;
  }[];
}

export function ModelGrid({ channel, models }: ModelGridProps) {
  const [query, setQuery] = useState("");
  const [visibleCount, setVisibleCount] = useState(PAGE_SIZE);
  const deferredQuery = useDeferredValue(query);

  const filtered = useMemo(() => {
    const raw = deferredQuery || "";
    const trimmed = raw.trim();
    if (!trimmed) return models;

    let normalizedQuery = trimmed.toLowerCase();
    if (normalizedQuery.endsWith(".mod")) {
      normalizedQuery = normalizedQuery.slice(0, -4);
    }

    if (normalizedQuery.length < 3) return models;

    return models.filter((model) => {
      const slug = model.slug.toLowerCase();
      const title = model.title.toLowerCase();
      return slug.includes(normalizedQuery) || title.includes(normalizedQuery);
    });
  }, [models, deferredQuery]);

  useEffect(() => {
    setVisibleCount(PAGE_SIZE);
  }, [deferredQuery]);

  const visible = filtered.slice(0, visibleCount);

  return (
    <div className="space-y-6">
      <div className="flex flex-col gap-3 sm:flex-row sm:items-center sm:justify-between">
        <label className="w-full sm:max-w-xs">
          <span className="text-xs font-semibold uppercase tracking-wide text-slate-500">Filter models</span>
          <input
            className="mt-1 w-full rounded-lg border border-slate-300 bg-white px-3 py-2 text-sm text-slate-900 transition focus:border-blue-500 focus:outline-none focus:ring-2 focus:ring-blue-500/50 dark:border-slate-700 dark:bg-slate-900 dark:text-slate-100"
            type="search"
            placeholder="Search by slug or title"
            value={query}
            onChange={(event) => setQuery(event.target.value)}
          />
        </label>
        <p className="text-xs text-slate-500 dark:text-slate-400">
          Showing <span className="font-semibold text-slate-700 dark:text-slate-200">{visible.length}</span> of
          <span className="font-semibold text-slate-700 dark:text-slate-200"> {filtered.length}</span> matching,
          <span className="font-semibold text-slate-700 dark:text-slate-200"> {models.length}</span> models
        </p>
      </div>

      <div className="grid gap-3 sm:grid-cols-2 lg:grid-cols-3">
        {visible.map((model) => (
          <Link
            key={model.slug}
            href={`/omnimodels/${encodeURIComponent(channel)}/${encodeURIComponent(model.slug)}`}
            prefetch={false}
            className="group rounded-lg border border-slate-200 bg-white p-4 text-sm shadow-sm transition hover:border-blue-400 hover:shadow-md dark:border-slate-800 dark:bg-slate-900"
          >
            <div className="flex items-center justify-between">
              <span className="font-semibold text-slate-800 transition group-hover:text-blue-600 dark:text-slate-100 group-hover:dark:text-blue-300">
                {model.title}
              </span>
              <span className="rounded-md bg-slate-100 px-2 py-1 text-xs font-mono text-slate-600 group-hover:bg-blue-500 group-hover:text-white dark:bg-slate-800 dark:text-slate-400">
                {model.slug}
              </span>
            </div>
          </Link>
        ))}
      </div>

      {visible.length < filtered.length && (
        <div className="flex justify-center">
          <button
            type="button"
            onClick={() => setVisibleCount((count) => Math.min(filtered.length, count + PAGE_SIZE))}
            className="rounded-lg border border-slate-300 bg-white px-5 py-2.5 text-sm font-medium text-slate-700 hover:bg-slate-50 dark:border-slate-700 dark:bg-slate-900 dark:text-slate-200 dark:hover:bg-slate-800"
          >
            Show {Math.min(PAGE_SIZE, filtered.length - visible.length)} more
          </button>
        </div>
      )}
    </div>
  );
}
